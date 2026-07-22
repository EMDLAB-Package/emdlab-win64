% EMDLAB: Electrical Machines Design Laboratory
% a two dimensional nonlinear magnetic-transient solver
% nonlinear solver: Newton-Raphson
% first order triangular mesh
% triangular lagrangian elements: 3 points per element
% isotropic
% homogenous
% nonlinear
% with motion: remesh technology with fixed number of nodes

classdef emdlab_solvers_mt2d_tl3_ihnlwm < handle & emdlab_solvers_mt2d_tlcp & matlab.mixin.Copyable

    properties
        % this flag is used to detect any motion in problem, when we have motion
        % we have to rebuild matrices
        isNeededToRebuild = false;
    end

    methods

        % initialization
        function obj = emdlab_solvers_mt2d_tl3_ihnlwm(m)

            % mesh pointer
            m.ggmesh(true);
            obj.m = m;

            % default settings for solver
            obj.solverSettings.relativeError = 1e-8;
            obj.solverSettings.maxIteration = 100;
            obj.solverSettings.relativeEnergyResidual = 1e-3;

            % set default properties of mesh zones
            mzNames = fieldnames(obj.m.mzs);

            for i = 1:numel(mzNames)
                obj.setdp(mzNames{i});
            end

        end

        % solver settings
        function setSolverMaxIteration(obj, maxIteration)

            if maxIteration < 0 || rem(maxIteration, 1)
                error('maxIteration must be a positive integer.');
            end

            obj.solverSettings.maxIteration = maxIteration;

        end

        function setSolverRelativeError(obj, relativeError)

            if relativeError < 0
                error('relativeError must be a real positive number.');
            end

            obj.solverSettings.relativeError = relativeError;

        end

        function setSolverRelativeEnergyResidual(obj, relativeEnergyResidual)

            if relativeEnergyResidual < 0
                error('relativeEnergyResidual must be a real positive number.');
            end

            obj.solverSettings.relativeEnergyResidual = relativeEnergyResidual;

        end

        function setMonitor(obj, value)

            % value = true or false
            obj.monitorResiduals = value;

        end

        % assign elements data
        function assignEdata(obj, initialRelativePermeability)

            % check states
            if obj.isElementDataAssigned, return; end

            % preparing mesh data
            obj.m.evalKeMeFe_TL3(true);
            tic, disp('-------------------------------------------------------');

            % assigning material and force data to each triangle for each mesh zone
            % getting mesh zones
            mzsName = fieldnames(obj.m.mzs);

            obj.edata.areAllLinear = true;
            % loop over mesh zones
            for i = 1:obj.m.Nmzs

                mzptr = obj.m.mzs.(mzsName{i});

                % allocate memory
                mzptr.props.MagneticReluctivity = zeros(1,mzptr.Ne);
                mzptr.props.ElectricConductivity = zeros(1,mzptr.Ne);
                mzptr.props.InternalCurrentDensity = zeros(1,mzptr.Ne);
                mzptr.props.MagnetizationX = zeros(1,mzptr.Ne);
                mzptr.props.MagnetizationY = zeros(1,mzptr.Ne);

                if ~obj.m.mts.(mzptr.material).MagneticPermeability.isIsotropic

                    throw(MException('', 'Some materials are non-isotropic, please select the correct solver.'));

                elseif obj.m.mts.(mzptr.material).MagneticPermeability.isLinear

                    % assigning Magnetic Permeability
                    mzptr.props.MagneticReluctivity(:) = 1/obj.m.mts.(mzptr.material).MagneticPermeability.value;

                else

                    obj.edata.areAllLinear = false;
                    if nargin == 2
                        mzptr.props.MagneticReluctivity(:) = obj.pcts.nu0 / initialRelativePermeability;
                    else
                        mzptr.props.MagneticReluctivity(:) = obj.pcts.nu0 / 500;
                    end

                end

                % assigning Electric Conductivity for activated zones for eddy currents
                if mzptr.props.isEddyZone
                    mzptr.props.ElectricConductivity(:) = obj.m.mts.(mzptr.material).ElectricConductivity.value;
                end

                % assigning Magnetization
                if mzptr.props.isMagnetized

                    M = mzptr.props.magnetization.getM(mzptr.getCenterOfElements);
                    mzptr.props.MagnetizationX(:) = M(:, 1)';
                    mzptr.props.MagnetizationY(:) = M(:, 2)';

                end

            end

            % evaluation of P and Q matrices of coils
            coilNames = fieldnames(obj.coils);

            for i = 1:obj.Ncoils
                % get coil pointer
                cptr = obj.coils.(coilNames{i});

                % initialize coil Pstranded, Psolid, Qvec, and R matrices for each coil arm
                cptr.Qvec = sparse(1,obj.m.Nn);
                for j = 1:cptr.NcoilArms
                    % pointer to coil arm
                    mzptr = obj.m.mzs.(cptr.coilArms(j));

                    % coil arm turns density
                    k = mzptr.props.turns/mzptr.getArea;

                    % for strandeds
                    % initialize value matrix
                    val = zeros(3,obj.m.Ne);
                    val(:,obj.m.ezi(:,mzptr.zi)) = k * obj.m.mtcs.Fe(:,obj.m.ezi(:,mzptr.zi));

                    % store matrices
                    mzptr.props.Pstranded = sparse(obj.m.cl', ones(3, obj.m.Ne), val);
                    mzptr.props.Qvec = mzptr.props.Pstranded' * obj.getDepth;
                    cptr.Qvec = cptr.Qvec + mzptr.props.direction * mzptr.props.Qvec;

                    % for solids
                    % initialize value matrix
                    val = zeros(3,obj.m.Ne);
                    val(:,obj.m.ezi(:,mzptr.zi)) = (obj.m.mts.(mzptr.material).ElectricConductivity.value * obj.units.k_length^2/obj.getDepth) * ...
                        obj.m.mtcs.Fe(:,obj.m.ezi(:,mzptr.zi));
                    mzptr.props.Psolid = sparse(obj.m.cl', ones(3,obj.m.Ne), val);

                end


            end

            disp('Initialization of material and force data compeleted.')
            toc, disp('-------------------------------------------------------');

            % Construction coils related matrices
            tic, disp('-------------------------------------------------------');

            % check windings
            obj.checkCoils;

            % allocating memory for coil related matrices
            obj.mtcs.K21 = zeros(obj.NcoilArms,obj.m.Nn);
            obj.mtcs.K31 = zeros(obj.Ncoils,obj.m.Nn);
            obj.mtcs.K12 = zeros(obj.m.Nn,obj.NcoilArms);
            obj.mtcs.K22 = zeros(obj.NcoilArms,obj.NcoilArms);
            obj.mtcs.K32 = zeros(obj.Ncoils,obj.NcoilArms);
            obj.mtcs.K13 = zeros(obj.m.Nn,obj.Ncoils);
            obj.mtcs.K23 = zeros(obj.NcoilArms,obj.Ncoils);
            obj.mtcs.K33 = zeros(obj.Ncoils,obj.Ncoils);
            obj.mtcs.Ksx = zeros(obj.NstarConnections,obj.m.Nn+obj.NcoilArms+obj.Ncoils);
            obj.mtcs.Ksy = zeros(obj.m.Nn+obj.NcoilArms+obj.Ncoils,obj.NstarConnections);
            obj.mtcs.Kss = zeros(obj.NstarConnections);

            % construction of matrices
            for i = 1:obj.Ncoils
                % get coil pointer
                cptr = obj.coils.(coilNames{i});
                % set coil current-voltage equation; the unknown is coil current
                if ~cptr.isCageMember
                    if cptr.isCurrentFed
                        % force current fed
                        obj.mtcs.K33(cptr.ci, cptr.ci) = 1;
                    else
                        % set coil Qvec to calculate total coil flux linkage
                        obj.mtcs.K31(cptr.ci,:) = cptr.Qvec;

                        % set coil resistance
                        obj.mtcs.K33(cptr.ci, cptr.ci) = cptr.Rdc;
                    end
                end
                % adjust matrices: importand K12 matrix
                if cptr.isStranded
                    % set Pvector of stranded mesh zones
                    for j = 1:cptr.NcoilArms
                        % pointer to coil arm
                        mzptr = obj.m.mzs.(cptr.coilArms{j});
                        obj.mtcs.K12(:,mzptr.props.cai) = -mzptr.props.Pstranded;
                        obj.mtcs.K22(mzptr.props.cai,mzptr.props.cai) = 1;
                        obj.mtcs.K23(mzptr.props.cai,cptr.ci) =  -mzptr.props.direction;
                    end
                else
                    % set Pvector of solid mesh zones
                    for j = 1:cptr.NcoilArms
                        % pointer to coil arm
                        mzptr = obj.m.mzs.(cptr.coilArms{j});
                        obj.mtcs.K12(:,mzptr.props.cai) = -mzptr.props.Psolid;
                        obj.mtcs.K21(mzptr.props.cai,:) = mzptr.props.Qvec;
                        obj.mtcs.K22(mzptr.props.cai,mzptr.props.cai) = -1;
                        obj.mtcs.K23(mzptr.props.cai,cptr.ci) =  mzptr.props.direction * mzptr.props.Rdc;
                    end
                end
            end

            % set voltage-current equation of cages
            cageNames = fieldnames(obj.cages);
            for i = 1:obj.Ncages
                % get pointer to cage
                cptr = obj.cages.(cageNames{i});

                % set voltage coefficienct matrix
                cptr.updateKuKrKl;
                obj.mtcs.K32(cptr.ciStart:cptr.ciEnd, cptr.caiStart:cptr.caiEnd) = cptr.Ku;

            end

            % set voltage-current equation of star connections
            starConnectionNames = fieldnames(obj.starConnections);
            for i = 1:obj.NstarConnections
                % get pointer to star connection
                scptr = obj.starConnections.(starConnectionNames{i});

                % set voltage coefficienct matrix
                obj.mtcs.Ksy(scptr.ci + obj.m.Nn + obj.NcoilArms) = -1;
                obj.mtcs.Ksx(scptr.ci + obj.m.Nn + obj.NcoilArms) = 1;

            end

            disp('Coils related matrices are constructed.')
            toc, disp('-------------------------------------------------------');

            % initialize results with zero A
            obj.results.A = zeros(obj.m.Nn, 1);
            obj.results.VICoilArms = zeros(obj.NcoilArms, 1);
            obj.results.ICoils = zeros(obj.Ncoils, 1);

            % change states
            obj.isElementDataAssigned = true;

        end

        function rebuildKeMeFe(obj)

            % Rebuild of [K], [M] and [Fm]
            tic, disp('-------------------------------------------------------');

            mzNames = fieldnames(obj.m.mzs);
            xNe = 0;
            for mzName = obj.m.getMeshZoneNames
                xNe = xNe + obj.m.mzs.(mzName).Ne;
            end

            obj.m.cl = zeros(xNe,3);
            obj.m.JIT = zeros(4,xNe);
            obj.m.mtcs.Ke = zeros(9,xNe);
            obj.m.mtcs.Me = zeros(9,xNe);
            obj.m.mtcs.Fe = zeros(3,xNe);
            obj.m.mtcs.FeMx = zeros(3,xNe);
            obj.m.mtcs.FeMy = zeros(3,xNe);

            obj.edata.MagneticReluctivity = zeros(1,xNe);
            obj.edata.ElectricConductivity = zeros(1,xNe);
            obj.edata.InternalCurrentDensity = zeros(1,xNe);
            obj.edata.MagnetizationX = zeros(1,xNe);
            obj.edata.MagnetizationY = zeros(1,xNe);

            index1 = 1;
            for i = 1:obj.m.Nmzs

                mzptr = obj.m.mzs.(mzNames{i});
                index2 = index1 + mzptr.Ne - 1;

                obj.m.cl(index1:index2,:) = mzptr.props.cl;
                obj.m.JIT(:,index1:index2) = mzptr.props.JIT;
                obj.m.mtcs.Ke(:,index1:index2) = mzptr.props.Ke;
                obj.m.mtcs.Me(:,index1:index2) = mzptr.props.Me;
                obj.m.mtcs.Fe(:,index1:index2) = mzptr.props.Fe;
                obj.m.mtcs.FeMx(:,index1:index2) = mzptr.props.FeMx;
                obj.m.mtcs.FeMy(:,index1:index2) = mzptr.props.FeMy;

                obj.edata.MagneticReluctivity(:,index1:index2) = mzptr.props.MagneticReluctivity;
                obj.edata.ElectricConductivity(:,index1:index2) = mzptr.props.ElectricConductivity;
                obj.edata.InternalCurrentDensity(:,index1:index2) = mzptr.props.InternalCurrentDensity;
                obj.edata.MagnetizationX(:,index1:index2) = mzptr.props.MagnetizationX;
                obj.edata.MagnetizationY(:,index1:index2) = mzptr.props.MagnetizationY;

                index1 = index2 + 1;

            end

            obj.m.elements = zeros(obj.m.Ne, 4);
            eindex = 0;

            % loop over mesh zones: for insertion of nodes and elements
            for i = 1:obj.m.Nmzs

                % specefying zone index
                obj.m.elements(1 + eindex:obj.m.mzs.(mzNames{i}).Ne + eindex, 4) = obj.m.mzs.(mzNames{i}).zi;
                eindex = eindex + obj.m.mzs.(mzNames{i}).Ne;

            end
            obj.m.setdata;
            obj.m.evalezi;

            % Assembeling [Fm]
            % assembling the load vector due to magnets
            obj.mtcs.Fm = (obj.edata.MagnetizationX .* obj.m.mtcs.FeMx + obj.edata.MagnetizationY .* obj.m.mtcs.FeMy) * (obj.units.k_length * obj.units.k_magnetisation);
            obj.mtcs.Fm = sparse(obj.m.cl', ones(3 * obj.m.Ne, 1), obj.mtcs.Fm);

            % Assembeling [K]
            [Iindex, Jindex] = emdlab_flib_getij(3,1);
            Iindex = obj.m.cl(:, Iindex)';
            Jindex = obj.m.cl(:, Jindex)';
            obj.mtcs.K11 = sparse(Iindex, Jindex, obj.edata.MagneticReluctivity .* obj.m.mtcs.Ke);
            obj.mtcs.M11 = sparse(Iindex, Jindex, obj.edata.ElectricConductivity .* obj.m.mtcs.Me * obj.units.k_length^2);

            disp('Rebuild of [K], [M], and [Fm] compeleted.');
            toc, disp('-------------------------------------------------------');

            obj.isNeededToRebuild = false;

        end

        function solveForInitialConditions(obj)

            % check if it is already solved for initial conditions
            if obj.isSolvedForInitialConditions, return; end

            % prerequisties
            obj.assignEdata;

            % updating boundary conditions
            obj.bcs.updateAll;

            % rebuild matrices in the case motion
            obj.rebuildKeMeFe;

            % Assembeling [F]
            F1 = obj.mtcs.Fm;
            F2 = zeros(obj.NcoilArms,1);
            F3 = zeros(obj.Ncoils,1);

            % construct [F3] and [F2]
            coilNames = fieldnames(obj.coils);
            K21 = obj.mtcs.K21;
            K31 = obj.mtcs.K31;
            K33 = obj.mtcs.K33;

            % adjust coil related matrices
            for i = 1:obj.Ncoils
                % get coil pointer
                cptr = obj.coils.(coilNames{i});
                % set coil current equation
                if cptr.isCurrentFed
                    % force stranded coil current
                    F3(cptr.ci) = cptr.getCurrent(obj.simTime(end));
                else
                    % set the arm coil voltage in such a way that we get the initial current
                    for j = 1:cptr.NcoilArms
                        % pointer to coil arm
                        mzptr = obj.m.mzs.(cptr.coilArms{j});
                        K21(mzptr.props.cai,:) = 0;
                    end
                    % force stranded coil voltage
                    K31(cptr.ci,:) = 0;
                    K33(cptr.ci,cptr.ci) = 1;
                    F3(cptr.ci) = cptr.initialCurrent;
                end
            end

            % construction of field circuit equations
            K = [obj.mtcs.K11, obj.mtcs.K12, obj.mtcs.K13
                K21, obj.mtcs.K22, obj.mtcs.K23
                K31, obj.mtcs.K32, K33];
            F = [F1;F2;F3];

            % imposing boundary conditions on [K] and [F]
            % dbcs
            if obj.bcs.Nd
                F(obj.bcs.iD) = obj.bcs.vD;
                K(obj.bcs.iD, :) = sparse(1:obj.bcs.Ndbcs, obj.bcs.iD, ones(1, obj.bcs.Ndbcs), obj.bcs.Ndbcs, obj.m.Nn+obj.NcoilArms+obj.Ncoils);
            end

            % opbcs
            if obj.bcs.Nop
                F(obj.bcs.mOP) = F(obj.bcs.mOP) - F(obj.bcs.sOP);
                F(obj.bcs.sOP) = 0;
                K(obj.bcs.mOP, :) = K(obj.bcs.mOP, :) - K(obj.bcs.sOP, :);
                K(obj.bcs.sOP, :) = sparse([1:obj.bcs.Nopbcs, 1:obj.bcs.Nopbcs], ...
                    [obj.bcs.mOP; obj.bcs.sOP], ones(1, 2 * obj.bcs.Nopbcs), obj.bcs.Nopbcs, obj.m.Nn+obj.NcoilArms+obj.Ncoils);
            end

            % epbcs
            if obj.bcs.Nep
                F(obj.bcs.mEP) = F(obj.bcs.mEP) + F(obj.bcs.sEP);
                F(obj.bcs.sEP) = 0;
                K(obj.bcs.mEP, :) = K(obj.bcs.mEP, :) + K(obj.bcs.sEP, :);
                K(obj.bcs.sEP, :) = sparse([1:obj.bcs.Nepbcs, 1:obj.bcs.Nepbcs], ...
                    [obj.bcs.mEP; obj.bcs.sEP], [ones(1, obj.bcs.Nepbcs), -ones(1, obj.bcs.Nepbcs)], obj.bcs.Nepbcs, obj.m.Nn+obj.NcoilArms+obj.Ncoils);
            end

            disp('All boundary condition imposed.');
            toc, disp('-------------------------------------------------------');

            % solving [K][U] = [F]
            tic, disp('-------------------------------------------------------');

            solVector = full(K\F);
            obj.results.A = solVector(1:obj.m.Nn);
            obj.results.VICoilArms = solVector(obj.m.Nn+1:obj.m.Nn+obj.NcoilArms);
            obj.results.ICoils = solVector(obj.m.Nn+obj.NcoilArms+1:obj.m.Nn+obj.NcoilArms+obj.Ncoils);

            for i = 1:obj.Ncoils

                % get coil pointer
                cptr = obj.coils.(coilNames{i});
                cptr.current(end) = solVector(obj.m.Nn + obj.NcoilArms + cptr.ci);

                % calculate and store coil flux linkage
                cptr.fluxLinkage(end) = cptr.Qvec * obj.results.A;

            end

            % calculate flux density in all elements
            obj.evalBe;

            % getting mesh zones
            mzsName = string(fieldnames(obj.m.mzs)');

            % save field data for core loss & magnet eddy current loss calculation
            for mzName = mzsName

                mzptr = obj.m.mzs.(mzName);
                if mzptr.props.isCoreLossActivated
                    mzptr.props.Bxg = obj.results.Bxg(obj.m.ezi(:,mzptr.zi));
                    mzptr.props.Byg = obj.results.Byg(obj.m.ezi(:,mzptr.zi));
                end

                if mzptr.props.isEddyZone
                    mzptr.props.Az = obj.results.A(mzptr.l2g);
                end

            end

            disp('initial geuss evaluated.')
            toc, disp('-------------------------------------------------------');

            % loop for nonlinear solver
            tic, disp('-------------------------------------------------------');

            % initials values
            EResidual = inf;
            RelError = inf;
            Iterations = 0;
            xNgt = obj.m.Ne;
            xNgp = obj.m.Nn+obj.NcoilArms+obj.Ncoils;

            [Iindex, Jindex] = emdlab_flib_getij(3,1);
            Iindex = obj.m.cl(:, Iindex)';
            Jindex = obj.m.cl(:, Jindex)';

            % preparing error monitoring
            if obj.monitorResiduals

                ERF = gcf; cla;
                set(ERF, 'Name', 'mt2d_tl3_ihnlwm solver', 'WindowStyle', 'Normal');
                er = animatedline('color', 'r', 'Linewidth', 1.2, 'Marker', 's', 'MarkerEdgeColor','k');
                title("Progress: " + num2str(0) + "%");
                ylabel('log10(||dA||/||A||)');
                xlabel('Iteration Number');
                cAxis = gca;
                box on;
                grid on;
                set(gca, 'ylim', [log10(obj.solverSettings.relativeError)-1,0]);
                set(gca, 'box', 'on');
                cAxis.XAxis.FontSize = 12;
                cAxis.YAxis.FontSize = 12;
                cAxis.Title.FontSize = 12;
                cAxis.YTick = log10(obj.solverSettings.relativeError)-1:1;
                cAxis.XLim(1) = 0;
                cAxis.XMinorGrid = 'on';
                cAxis.YMinorGrid = 'off';

            end

            % solver history
            obj.solverHistory.relativeError = [];
            obj.solverHistory.totalEnergy = [];
            obj.solverHistory.totalConergy = [];

            % memory allocation for dnudB2
            dnudB2 = zeros(1, xNgt);

            % inintial value of alphaNR
            alphaNR = 0.7;

            % loop for non-linearity
            fprintf('Iter|Error   |Residual|time\n');
            while ((RelError > obj.solverSettings.relativeError) || (EResidual>obj.solverSettings.relativeEnergyResidual)) && ...
                    (Iterations < obj.solverSettings.maxIteration) && ~obj.edata.areAllLinear && any(F)

                % starting loop time
                loopTime = tic;

                % evaluation of B2 for each elements
                obj.evalBe;
                [obj.solverHistory.totalEnergy(end + 1),obj.solverHistory.totalConergy(end + 1)] = obj.evalTotalEnergyCoenergy;
                Bk = obj.results.Bxg.^2 + obj.results.Byg.^2;

                % calc energy residual
                if length(obj.solverHistory.totalEnergy)>2
                    EResidual = abs(obj.solverHistory.totalEnergy(end)-obj.solverHistory.totalEnergy(end-1))/...
                        obj.solverHistory.totalEnergy(end);
                end

                % updating nu & dnudB2
                for i = 1:obj.m.Nmts
                    mtptr = obj.m.mts.(obj.m.materialNames(i));

                    if ~mtptr.MagneticPermeability.isLinear
                        obj.edata.MagneticReluctivity(obj.m.emi(i,:)) = ppval(mtptr.vB2, Bk(obj.m.emi(i,:)));
                        dnudB2(obj.m.emi(i,:)) = ppval(mtptr.dvdB2, Bk(obj.m.emi(i,:)));
                    end

                end

                % construction of stiffness matrix [K]
                K11 = sparse(Iindex, Jindex, obj.edata.MagneticReluctivity .* obj.m.mtcs.Ke);

                K = [K11, obj.mtcs.K12, obj.mtcs.K13
                    K21, obj.mtcs.K22, obj.mtcs.K23
                    K31, obj.mtcs.K32, K33];

                % construction of [K] and [F] in NR algorithm
                FF = -K * solVector + F;

                % evaluation and adding of jacobian matrix
                K11 = K11 + sparse(Iindex, Jindex, emdlab_m2d_tl3_evalG(obj.m.cl, obj.m.mtcs.Ke, obj.m.JIT, obj.results.A, dnudB2) / obj.units.k_length^2);

                K = [K11, obj.mtcs.K12, obj.mtcs.K13
                    K21, obj.mtcs.K22, obj.mtcs.K23
                    K31, obj.mtcs.K32, K33];

                % imposing boundary conditions on incrimentals
                % dbcs
                if obj.bcs.Nd
                    FF(obj.bcs.iD) = obj.bcs.vD;
                    K(obj.bcs.iD, :) = sparse(1:obj.bcs.Ndbcs, obj.bcs.iD, ones(1, obj.bcs.Ndbcs), obj.bcs.Ndbcs, xNgp);
                end

                % opbcs
                if obj.bcs.Nop
                    FF(obj.bcs.mOP) = FF(obj.bcs.mOP) - FF(obj.bcs.sOP);
                    FF(obj.bcs.sOP) = 0;
                    K(obj.bcs.mOP, :) = K(obj.bcs.mOP, :) - K(obj.bcs.sOP, :);
                    K(obj.bcs.sOP, :) = sparse([1:obj.bcs.Nopbcs, 1:obj.bcs.Nopbcs], ...
                        [obj.bcs.mOP; obj.bcs.sOP], ones(1, 2 * obj.bcs.Nopbcs), obj.bcs.Nopbcs, xNgp);
                end

                % epbcs
                if obj.bcs.Nep
                    FF(obj.bcs.mEP) = FF(obj.bcs.mEP) + FF(obj.bcs.sEP);
                    FF(obj.bcs.sEP) = 0;
                    K(obj.bcs.mEP, :) = K(obj.bcs.mEP, :) + K(obj.bcs.sEP, :);
                    K(obj.bcs.sEP, :) = sparse([1:obj.bcs.Nepbcs, 1:obj.bcs.Nepbcs], ...
                        [obj.bcs.mEP; obj.bcs.sEP], [ones(1, obj.bcs.Nepbcs), -ones(1, obj.bcs.Nepbcs)], obj.bcs.Nepbcs, xNgp);
                end

                % solving [K][U] = [F]
                dU = full(K\FF);
                solVector = solVector + alphaNR*dU;

                obj.results.A = solVector(1:obj.m.Nn);
                obj.results.VICoilArms = solVector(obj.m.Nn+1:obj.m.Nn+obj.NcoilArms);
                obj.results.ICoils = solVector(obj.m.Nn+obj.NcoilArms+1:obj.m.Nn+obj.NcoilArms+obj.Ncoils);

                for i = 1:obj.Ncoils
                    % get coil pointer
                    cptr = obj.coils.(coilNames{i});
                    cptr.current(end) = solVector(obj.m.Nn + obj.NcoilArms + cptr.ci);
                    % calculate and store coil flux linkage
                    cptr.fluxLinkage(end) = cptr.Qvec * obj.results.A;
                end

                % check for convergency
                Residual = norm(dU, 2);
                RelError = Residual / norm(solVector, 2);

                % monitoring of error
                if obj.monitorResiduals
                    addpoints(er, Iterations+1, log10(RelError));
                    cAxis.XLim(2) = Iterations+2;
                    drawnow;
                end

                % solver history
                obj.solverHistory.relativeError(end + 1) = RelError;

                % printing Residual and RelError
                fprintf('->%2d|%.2e|%.2e|%0.3f\n', Iterations, RelError, Residual, toc(loopTime));

                % go to next iteration
                Iterations = Iterations + 1;

                % update alphaNR
                if length(obj.solverHistory.relativeError)>2
                    if obj.solverHistory.relativeError(end) > 0.8*obj.solverHistory.relativeError(end-1)
                        alphaNR = max((0.95-2e-2*rand)*alphaNR,0.5);
                    else
                        alphaNR = min((1.05+2e-2*rand)*alphaNR,0.9);
                    end
                end

            end

            if obj.monitorResiduals
                cAxis.YLim(2) = ceil(log10(obj.solverHistory.relativeError(1)));
            end

            obj.solverHistory.iterations = Iterations;
            disp(['Number of total iterations = ', num2str(Iterations - 1)]);
            toc, disp('-------------------------------------------------------');

            % update field quantities
            obj.evalBe;
            obj.evalHe;
            obj.evalBn;
            obj.evalHn;

            % save field data for core loss calculation & update magnetic reluctivities
            for mzName = mzsName

                mzptr = obj.m.mzs.(char(mzName));
                if mzptr.props.isCoreLossActivated
                    mzptr.props.Bxg(end,:) = obj.results.Bxg(obj.m.ezi(:,mzptr.zi));
                    mzptr.props.Byg(end,:) = obj.results.Byg(obj.m.ezi(:,mzptr.zi));
                end

                if mzptr.props.isEddyZone
                    mzptr.props.Az(:,end) = obj.results.A(mzptr.l2g);
                end

                mzptr.props.MagneticReluctivity = obj.edata.MagneticReluctivity(obj.m.ezi(:,mzptr.zi));

            end

            % change states
            obj.isSolvedForInitialConditions = true;

        end

        % solver core
        function obj = solveForOneTimeStep(obj, DeltaTime)

            % prerequisties
            obj.solveForInitialConditions;
            obj.simTime(end+1) = obj.simTime(end) + DeltaTime;
            obj.assignEdata;

            % updating boundary conditions
            obj.bcs.updateAll;

            if obj.isNeededToRebuild
                obj.rebuildKeMeFe;
            end

            % Assembeling [F]
            F1 = obj.mtcs.Fm + obj.mtcs.M11 * obj.results.A / DeltaTime;
            F2 = zeros(obj.NcoilArms,1);
            F3 = zeros(obj.Ncoils,1);
            F4 = zeros(obj.NstarConnections,1);

            % construct [F3] and [F2]
            coilNames = fieldnames(obj.coils);
            K31 = obj.mtcs.K31;
            K21 = obj.mtcs.K21;
            K33 = obj.mtcs.K33;

            for i = 1:obj.Ncoils

                % get coil pointer
                cptr = obj.coils.(coilNames{i});

                % adjust matrices
                % set coil current equation
                if ~cptr.isCageMember
                    switch cptr.fedType
                        case 'current'

                            % force stranded coil current
                            F3(cptr.ci) = cptr.getCurrent(obj.simTime(end));

                        case 'voltage'

                            % force stranded coil voltage
                            K31(cptr.ci,:) = K31(cptr.ci,:)/DeltaTime;
                            F3(cptr.ci) = cptr.getVoltage(obj.simTime(end)) + K31(cptr.ci,:) * obj.results.A;
                    end
                end

                if strcmpi(cptr.eddyType, 'solid')
                    for j = 1:cptr.NcoilArms

                        % pointer to coil arm
                        mzptr = obj.m.mzs.(cptr.coilArms(j));

                        K21(mzptr.props.cai,:) = K21(mzptr.props.cai,:)/DeltaTime;

                        F2(mzptr.props.cai) =  K21(mzptr.props.cai,:) * obj.results.A;

                    end
                end

            end

            % set voltage-current equation of cages
            cageNames = fieldnames(obj.cages);
            for i = 1:obj.Ncages
                % get pointer to cage
                cptr = obj.cages.(cageNames{i});

                % set current coefficienct matrix
                K33(cptr.ciStart:cptr.ciEnd, cptr.ciStart:cptr.ciEnd) = cptr.Kr + cptr.Kl / DeltaTime;
                F3(cptr.ciStart:cptr.ciEnd) = cptr.Kl * obj.results.ICoils(cptr.ciStart:cptr.ciEnd) / DeltaTime;
            end

            % construction of field circuit equations
            K = [obj.mtcs.K11 + obj.mtcs.M11/DeltaTime, obj.mtcs.K12, obj.mtcs.K13
                K21, obj.mtcs.K22, obj.mtcs.K23
                K31, obj.mtcs.K32, K33];
            K = [K,obj.mtcs.Ksy
                obj.mtcs.Ksx,obj.mtcs.Kss];
            F = [F1;F2;F3;F4];

            tic, disp('-------------------------------------------------------');

            % imposing boundary conditions on [K] and [F]
            % dbcs
            if obj.bcs.Nd
                F(obj.bcs.iD) = obj.bcs.vD;
                K(obj.bcs.iD, :) = sparse(1:obj.bcs.Ndbcs, obj.bcs.iD, ones(1, obj.bcs.Ndbcs), obj.bcs.Ndbcs, obj.m.Nn+obj.NcoilArms+obj.Ncoils+obj.NstarConnections);
            end

            % opbcs
            if obj.bcs.Nop
                F(obj.bcs.mOP) = F(obj.bcs.mOP) - F(obj.bcs.sOP);
                F(obj.bcs.sOP) = 0;
                K(obj.bcs.mOP, :) = K(obj.bcs.mOP, :) - K(obj.bcs.sOP, :);
                K(obj.bcs.sOP, :) = sparse([1:obj.bcs.Nopbcs, 1:obj.bcs.Nopbcs], ...
                    [obj.bcs.mOP; obj.bcs.sOP], ones(1, 2 * obj.bcs.Nopbcs), obj.bcs.Nopbcs, obj.m.Nn+obj.NcoilArms+obj.Ncoils+obj.NstarConnections);
            end

            % epbcs
            if obj.bcs.Nep
                F(obj.bcs.mEP) = F(obj.bcs.mEP) + F(obj.bcs.sEP);
                F(obj.bcs.sEP) = 0;
                K(obj.bcs.mEP, :) = K(obj.bcs.mEP, :) + K(obj.bcs.sEP, :);
                K(obj.bcs.sEP, :) = sparse([1:obj.bcs.Nepbcs, 1:obj.bcs.Nepbcs], ...
                    [obj.bcs.mEP; obj.bcs.sEP], [ones(1, obj.bcs.Nepbcs), -ones(1, obj.bcs.Nepbcs)], obj.bcs.Nepbcs, obj.m.Nn+obj.NcoilArms+obj.Ncoils+obj.NstarConnections);
            end

            disp('All boundary condition imposed.');
            toc, disp('-------------------------------------------------------');

            % solving [K][U] = [F]
            tic, disp('-------------------------------------------------------');

            solVector = full(K \ F);
            obj.results.A = solVector(1:obj.m.Nn);
            obj.results.VICoilArms = solVector(obj.m.Nn+1:obj.m.Nn+obj.NcoilArms);
            obj.results.ICoils = solVector(obj.m.Nn+obj.NcoilArms+1:obj.m.Nn+obj.NcoilArms+obj.Ncoils);

            for i = 1:obj.Ncoils
                % get coil pointer
                cptr = obj.coils.(coilNames{i});
                cptr.current(end+1) = solVector(obj.m.Nn + obj.NcoilArms + cptr.ci);
                cptr.fluxLinkage(end+1) = cptr.Qvec * obj.results.A;
                cptr.inducedVoltage(end+1) = (cptr.fluxLinkage(end)-cptr.fluxLinkage(end-1))/DeltaTime;
                %                 if length(cptr.fluxLinkage)>2
                %                     cptr.inducedVoltage(end+1) = (cptr.fluxLinkage(end)-cptr.fluxLinkage(end-2))/(2*DeltaTime);
                %                 else
                %                     cptr.inducedVoltage(end+1) = (cptr.fluxLinkage(end)-cptr.fluxLinkage(end-1))/DeltaTime;
                %                 end
                cptr.voltage(end+1) = cptr.Rdc * cptr.current(end) + cptr.inducedVoltage(end);
            end

            obj.evalBe;

            % getting mesh zones
            mzsName = string(fieldnames(obj.m.mzs)');

            % save field data for core loss calculation
            for mzName = mzsName

                mzptr = obj.m.mzs.(mzName);
                if mzptr.props.isCoreLossActivated
                    mzptr.props.Bxg(end+1,:) = obj.results.Bxg(obj.m.ezi(:,mzptr.zi));
                    mzptr.props.Byg(end+1,:) = obj.results.Byg(obj.m.ezi(:,mzptr.zi));
                end

                if mzptr.props.isEddyZone
                    mzptr.props.Az(:,end+1) = obj.results.A(mzptr.l2g);
                end

            end

            disp('initial geuss evaluated.')
            toc, disp('-------------------------------------------------------');

            % loop for nonlinear solver
            tic, disp('-------------------------------------------------------');

            % initials values
            RelEResidual = inf;
            RelError = inf;
            Iterations = 0;
            xNgt = obj.m.Ne;
            xNgp = obj.m.Nn+obj.NcoilArms+obj.Ncoils+obj.NstarConnections;

            [Iindex, Jindex] = emdlab_flib_getij(3,1);
            Iindex = obj.m.cl(:, Iindex)';
            Jindex = obj.m.cl(:, Jindex)';

            % preparing error monitoring
            if obj.monitorResiduals

                ERF = gcf; cla;
                set(ERF, 'Name', 'mt2d_tl3_ihnlwtm solver', 'WindowStyle', 'Normal');
                er = animatedline('color', 'r', 'Linewidth', 1.2, 'Marker', 's', 'MarkerEdgeColor','k');
                ylabel('log10(||dA||/||A||)');
                xlabel('Iteration Number');
                cAxis = gca;
                box on;
                grid on;
                title("Time: " + sprintf('%4.4f',obj.simTime(end)));
                set(gca, 'ylim', [log10(obj.solverSettings.relativeError)-1,0]);
                set(gca, 'box', 'on');
                cAxis.XAxis.FontSize = 12;
                cAxis.YAxis.FontSize = 12;
                cAxis.Title.FontSize = 12;
                cAxis.YTick = log10(obj.solverSettings.relativeError)-1:1;
                cAxis.XLim(1) = 0;
                cAxis.XMinorGrid = 'on';
                cAxis.YMinorGrid = 'off';

            end

            % solver history
            obj.solverHistory.relativeError = [];
            obj.solverHistory.totalEnergy = [];
            obj.solverHistory.totalConergy = [];

            % memory allocation for dnudB2
            dnudB2 = zeros(1, xNgt);

            % inintial value of alphaNR
            alphaNR = 0.7;

            % loop for non-linearity
            fprintf('Iter|Error   |Residual|time\n');
            while ((RelError > obj.solverSettings.relativeError) || (RelEResidual>obj.solverSettings.relativeEnergyResidual)) ...
                    && (Iterations < obj.solverSettings.maxIteration) && ~obj.edata.areAllLinear && any(F)

                % starting loop time
                loopTime = tic;

                % evaluation of B2 for each elements
                obj.evalBe;
                [obj.solverHistory.totalEnergy(end + 1),obj.solverHistory.totalConergy(end + 1)] = obj.evalTotalEnergyCoenergy;
                Bk = obj.results.Bxg.^2 + obj.results.Byg.^2;

                % calculate relative energy residual
                if length(obj.solverHistory.totalEnergy)>2
                    RelEResidual = abs(obj.solverHistory.totalEnergy(end)-obj.solverHistory.totalEnergy(end-1))/...
                        obj.solverHistory.totalEnergy(end);
                end

                % updating nu & dnudB2
                for i = 1:obj.m.Nmts
                    mtptr = obj.m.mts.(obj.m.materialNames(i));

                    if ~mtptr.MagneticPermeability.isLinear
                        obj.edata.MagneticReluctivity(obj.m.emi(i,:)) = ppval(mtptr.vB2, Bk(obj.m.emi(i,:)));
                        dnudB2(obj.m.emi(i,:)) = ppval(mtptr.dvdB2, Bk(obj.m.emi(i,:)));
                    end

                end

                % construction of stiffness matrix [K]
                K11 = sparse(Iindex, Jindex, obj.edata.MagneticReluctivity .* obj.m.mtcs.Ke);

                K = [K11 + obj.mtcs.M11/DeltaTime, obj.mtcs.K12, obj.mtcs.K13
                    K21, obj.mtcs.K22, obj.mtcs.K23
                    K31, obj.mtcs.K32, K33];
                K = [K,obj.mtcs.Ksy
                    obj.mtcs.Ksx,obj.mtcs.Kss];

                % construction of [K] and [F] in NR algorithm
                FF = -K * solVector + F;

                % evaluation and adding of jacobian matrix
                K11 = K11 + sparse(Iindex, Jindex, emdlab_m2d_tl3_evalG(obj.m.cl, obj.m.mtcs.Ke, obj.m.JIT, obj.results.A, dnudB2) / obj.units.k_length^2);

                K = [K11 + obj.mtcs.M11/DeltaTime, obj.mtcs.K12, obj.mtcs.K13
                    K21, obj.mtcs.K22, obj.mtcs.K23
                    K31, obj.mtcs.K32, K33];
                K = [K,obj.mtcs.Ksy
                    obj.mtcs.Ksx,obj.mtcs.Kss];

                % imposing boundary conditions on incrimentals
                % dbcs
                if obj.bcs.Nd
                    FF(obj.bcs.iD) = obj.bcs.vD;
                    K(obj.bcs.iD, :) = sparse(1:obj.bcs.Ndbcs, obj.bcs.iD, ones(1, obj.bcs.Ndbcs), obj.bcs.Ndbcs, xNgp);
                end

                % opbcs
                if obj.bcs.Nop
                    FF(obj.bcs.mOP) = FF(obj.bcs.mOP) - FF(obj.bcs.sOP);
                    FF(obj.bcs.sOP) = 0;
                    K(obj.bcs.mOP, :) = K(obj.bcs.mOP, :) - K(obj.bcs.sOP, :);
                    K(obj.bcs.sOP, :) = sparse([1:obj.bcs.Nopbcs, 1:obj.bcs.Nopbcs], ...
                        [obj.bcs.mOP; obj.bcs.sOP], ones(1, 2 * obj.bcs.Nopbcs), obj.bcs.Nopbcs, xNgp);
                end

                % epbcs
                if obj.bcs.Nep
                    FF(obj.bcs.mEP) = FF(obj.bcs.mEP) + FF(obj.bcs.sEP);
                    FF(obj.bcs.sEP) = 0;
                    K(obj.bcs.mEP, :) = K(obj.bcs.mEP, :) + K(obj.bcs.sEP, :);
                    K(obj.bcs.sEP, :) = sparse([1:obj.bcs.Nepbcs, 1:obj.bcs.Nepbcs], ...
                        [obj.bcs.mEP; obj.bcs.sEP], [ones(1, obj.bcs.Nepbcs), -ones(1, obj.bcs.Nepbcs)], obj.bcs.Nepbcs, xNgp);
                end

                % solving [K][U] = [F]
                dU = K \ FF;
                solVector = full(solVector + alphaNR*dU);
                obj.results.A = solVector(1:obj.m.Nn);

                % check for convergency
                Residual = norm(dU, 2);
                RelError = Residual / norm(solVector, 2);

                % monitoring of error
                if obj.monitorResiduals
                    addpoints(er, Iterations+1, log10(RelError));
                    cAxis.XLim(2) = Iterations+2;
                    drawnow;
                end

                % solver history
                obj.solverHistory.relativeError(end + 1) = RelError;

                % printing Residual and RelError
                fprintf('->%2d|%.2e|%.2e|%0.3f\n', Iterations, RelError, Residual, toc(loopTime));

                % go to next iteration
                Iterations = Iterations + 1;

                % update alphaNR
                if length(obj.solverHistory.relativeError)>2
                    if obj.solverHistory.relativeError(end) > 0.8*obj.solverHistory.relativeError(end-1)
                        alphaNR = max((0.95-2e-2*rand)*alphaNR,0.5);
                    else
                        alphaNR = min((1.05+2e-2*rand)*alphaNR,0.9);
                    end
                end

            end

            if obj.monitorResiduals
                cAxis.YLim(2) = ceil(log10(obj.solverHistory.relativeError(1)));
            end

            obj.solverHistory.iterations = Iterations;
            disp(['Number of total iterations = ', num2str(Iterations - 1)]);
            toc, disp('-------------------------------------------------------');

            % change states
            obj.evalBe;
            obj.evalHe;
            obj.evalBn;
            obj.evalHn;

            % save information of coils
            obj.results.VICoilArms = solVector(obj.m.Nn+1:obj.m.Nn+obj.NcoilArms);
            obj.results.ICoils = solVector(obj.m.Nn+obj.NcoilArms+1:obj.m.Nn+obj.NcoilArms+obj.Ncoils);

            for i = 1:obj.Ncoils
                % get coil pointer
                cptr = obj.coils.(coilNames{i});
                cptr.current(end) = solVector(obj.m.Nn + obj.NcoilArms + cptr.ci);
                cptr.fluxLinkage(end) = cptr.Qvec * obj.results.A;
                cptr.inducedVoltage(end) = (cptr.fluxLinkage(end)-cptr.fluxLinkage(end-1))/DeltaTime;
                %                     if length(cptr.fluxLinkage)>2
                %                         cptr.inducedVoltage(end) = (cptr.fluxLinkage(end)-cptr.fluxLinkage(end-2))/(2*DeltaTime);
                %                     else
                %                         cptr.inducedVoltage(end) = (cptr.fluxLinkage(end)-cptr.fluxLinkage(end-1))/DeltaTime;
                %                     end
                cptr.voltage(end) = cptr.Rdc * cptr.current(end) + cptr.inducedVoltage(end);
            end

            % save field data for core loss calculation & update magnetic reluctivities
            for mzName = mzsName

                mzptr = obj.m.mzs.(char(mzName));
                if mzptr.props.isCoreLossActivated
                    mzptr.props.Bxg(end,:) = obj.results.Bxg(obj.m.ezi(:,mzptr.zi));
                    mzptr.props.Byg(end,:) = obj.results.Byg(obj.m.ezi(:,mzptr.zi));
                end

                if mzptr.props.isEddyZone
                    mzptr.props.Az(:,end) = obj.results.A(mzptr.l2g);
                end

                mzptr.props.MagneticReluctivity = obj.edata.MagneticReluctivity(obj.m.ezi(:,mzptr.zi));

            end

        end

        function solve(obj, stopTime, timeStep)

            if stopTime <= obj.simTime(end)
                return;
            end

            timeInterval = stopTime - obj.simTime(end);
            spentTime = 0;
            DeltaTime = diff(obj.simTime(end):timeStep:stopTime);
            for dt = DeltaTime
                obj.solveForOneTimeStep(dt);
                spentTime = spentTime + dt;
                if obj.monitorResiduals
                    title("Progress: " + sprintf('%3.2f',100*spentTime/timeInterval) + "%");
                end
            end

        end

        function rotateMovingRegion(obj, movingRegionName, rotAngle, xc, yc)

            % moving region rotation
            tic, disp('-------------------------------------------------------');

            if nargin < 4
                xc = 0;
                yc = 0;
            end

            if ~obj.isSolvedForInitialConditions
                obj.assignEdata;
            end

            movingRegionName = obj.checkMovingRegionExistence(movingRegionName);

            % get pointer to moving region
            mrptr = obj.movingRegions.(movingRegionName);
            mrptr.motionHistory(end+1,:) = [0, 0, xc, yc, rotAngle];

            obj.rotateMeshZones(mrptr.meshZones, rotAngle, xc, yc);
            obj.m.imzs.(mrptr.interface).rotate(rotAngle, [xc, yc]);

            % set interface mesh zone new properties
            mzptr = obj.m.mzs.(mrptr.interface);
            mzptr.nodes = obj.m.imzs.(mrptr.interface).m.nodes;
            mzptr.cl = obj.m.imzs.(mrptr.interface).m.cl;
            mzptr.setdataForce;
            mzptr.props.cl = mzptr.l2g(mzptr.cl);
            [mzptr.props.JIT,mzptr.props.gea] = emdlab_m2d_tl3_evalJIT(mzptr.cl, mzptr.nodes);
            [mzptr.props.Ke, mzptr.props.Me, mzptr.props.Fe, mzptr.props.FeMx, mzptr.props.FeMy] = emdlab_m2d_tl3_evalKeMeFe(mzptr.cl, mzptr.nodes);
            mzptr.props.MagneticReluctivity = ones(1,mzptr.Ne)*obj.pcts.nu0;
            mzptr.props.ElectricConductivity = zeros(1,mzptr.Ne);
            mzptr.props.InternalCurrentDensity = zeros(1,mzptr.Ne);
            mzptr.props.MagnetizationX = zeros(1,mzptr.Ne);
            mzptr.props.MagnetizationY = zeros(1,mzptr.Ne);

            disp('rotation compeleted.');
            toc, disp('-------------------------------------------------------');

            % updateJIT for rotating mesh zones

            obj.isNeededToRebuild = true;

        end

    end
end
