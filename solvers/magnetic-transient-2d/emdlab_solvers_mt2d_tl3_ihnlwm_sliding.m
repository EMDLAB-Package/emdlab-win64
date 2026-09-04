% EMDLAB: Electrical Machines Design Laboratory
% two dimensional nonlinear magnetic-transient solver
% nonlinear solver: Newton-Raphson
% first order triangular mesh
% triangular lagrangian elements: 3 points per element
% isotropic
% homogenous
% nonlinear
% without motion

classdef emdlab_solvers_mt2d_tl3_ihnlwm_sliding < handle & emdlab_solvers_mt2d_tlcp

    methods

        function obj = emdlab_solvers_mt2d_tl3_ihnlwm_sliding(m)

            % generate global mesh & set mesh pointer
            m.ggmesh;
            obj.m = m;

            % default settings for solver
            obj.solverSettings.relativeError = 1e-8;
            obj.solverSettings.maxIteration = 100;
            obj.solverSettings.relativeEnergyResidual = 1e-3;

            % set default properties of mesh zones
            for mzName = obj.m.getMeshZoneNames
                obj.setdp(mzName);
            end

        end

        function assignEdata(obj, InitNur)
            % assign elements data
            % assigning material and force data to each triangle element

            % check states
            if obj.isElementDataAssigned, return; end

            % preparing mesh data
            obj.m.evalKeMeFe_TL3;
            timeHolder = tic;

            % allocation of memory
            obj.edata.MagneticReluctivity = zeros(1, obj.m.Ne);
            obj.edata.ElectricConductivity = zeros(1, obj.m.Ne);
            obj.edata.InternalCurrentDensity = zeros(1, obj.m.Ne);
            obj.edata.MagnetizationX = zeros(1, obj.m.Ne);
            obj.edata.MagnetizationY = zeros(1, obj.m.Ne);

            % getting mesh zones
            mzsName = fieldnames(obj.m.mzs);

            obj.edata.areAllLinear = true;
            % loop over mesh zones
            for i = 1:obj.m.Nmzs

                mzptr = obj.m.mzs.(mzsName{i});

                if ~obj.m.mts.(mzptr.material).MagneticPermeability.isIsotropic

                    throw(MException('', 'Some materials are non-isotropic, please select the correct solver.'));

                elseif obj.m.mts.(mzptr.material).MagneticPermeability.isLinear

                    % assigning Magnetic Permeability
                    obj.edata.MagneticReluctivity(obj.m.ezi(:, mzptr.zi)) = 1/obj.m.mts.(mzptr.material).MagneticPermeability.value;

                else

                    obj.edata.areAllLinear = false;
                    if nargin == 2
                        obj.edata.MagneticReluctivity(obj.m.ezi(:, mzptr.zi)) = InitNur * obj.pcts.nu0;
                    else
                        obj.edata.MagneticReluctivity(obj.m.ezi(:, mzptr.zi)) = 0.001 * obj.pcts.nu0;
                    end

                end

                % assigning Electric Conductivity for activated zones for eddy currents
                if mzptr.props.isEddyZone
                    obj.edata.ElectricConductivity(obj.m.ezi(:, mzptr.zi)) = obj.m.mts.(mzptr.material).ElectricConductivity.value;
                end

                % assigning Magnetization
                if mzptr.props.isMagnetized

                    M = mzptr.props.magnetization.getM(mzptr.getCenterOfElements);
                    obj.edata.MagnetizationX(obj.m.ezi(:, mzptr.zi)) = M(:, 1)';
                    obj.edata.MagnetizationY(obj.m.ezi(:, mzptr.zi)) = M(:, 2)';

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

            obj.dispMessageLine('Initialization of material and force data compeleted.', timeHolder);

            % Construction coils related matrices
            timeHolder = tic;

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

            obj.dispMessageLine('Coils related matrices are constructed.', timeHolder);

            % Construction of [K], [M] and [Fm]
            timeHolder = tic;

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
            obj.dispMessageLine('Construction of [K], [M], and [Fm] compeleted.', timeHolder);

            % initialize results with zero A
            obj.results.A = zeros(obj.m.Nn, 1);
            obj.results.VICoilArms = zeros(obj.NcoilArms, 1);
            obj.results.ICoils = zeros(obj.Ncoils, 1);

            % change states
            obj.isElementDataAssigned = true;

        end

        function solveForInitialConditions(obj)

            % check if it is already solved for initial conditions
            if obj.isSolvedForInitialConditions, return; end

            % prerequisties
            obj.assignEdata;

            % updating boundary conditions
            obj.bcs.updateAll;

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
            if ~any(F), return; end

            tic, disp('-------------------------------------------------------');

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

            % apply sliding contacts


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

            obj.evalBe;
            disp('initial geuss evaluated.')
            toc, disp('-------------------------------------------------------');

            if obj.edata.areAllLinear
                obj.evalHe;
                obj.evalBn;
                obj.evalHn;
                return
            end

            % loop for nonlinear solver
            tic, disp('-------------------------------------------------------');

            % initials values
            RelEResidual = inf;
            RelError = inf;
            Iterations = 0;
            xNgt = obj.m.Ne;
            xNgp = obj.m.Nn+obj.NcoilArms+obj.Ncoils;

            [Iindex, Jindex] = getij(3,1);
            Iindex = obj.m.cl(:, Iindex)';
            Jindex = obj.m.cl(:, Jindex)';

            % preparing error monitoring
            if obj.monitorResiduals

                ERF = gcf; cla;
                set(ERF, 'Name', 'mt2d_tl3_ihnlwtm solver', 'WindowStyle', 'Normal');
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
            while ((RelError > obj.solverSettings.relativeError) || (RelEResidual>obj.solverSettings.relativeEnergyResidual)) && (Iterations < obj.solverSettings.maxIteration)

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

                % apply sliding contacts


                % solving [K][U] = [F]
                dU = full(K\FF);
                solVector = solVector + alphaNR*dU;

                obj.results.A = solVector(1:obj.m.Nn);
                obj.results.VICoilArms = solVector(obj.m.Nn+1:obj.m.Nn+obj.NcoilArms);

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
                        alphaNR = min((1.05+2e-2*rand)*alphaNR+2e-3*rand,0.9);
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

            % change states
            obj.isSolvedForInitialConditions = true;

        end

        function obj = solveForOneTimeStep(obj, DeltaTime)

            obj.solveForInitialConditions;
            obj.simTime(end+1) = obj.simTime(end) + DeltaTime;

            % prerequisties
            obj.assignEdata;

            % updating boundary conditions
            obj.bcs.updateAll;

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

            % apply sliding contacts
            mi = obj.m.contacts.ag.m;
            si = obj.m.contacts.ag.s;

            a1 = atan_02pi(obj.m.nodes(mi,:));
            a2 = atan_02pi(obj.m.nodes(si,:));

            [a1,idx] = sort(a1);
            mi = mi(idx);
            [a2,idx] = sort(a2);
            si = si(idx);

            for i = 1:length(mi)

                j = 1;
                while (a1(i) > a2(j)) && (j <= length(a2))
                    j = j + 1;
                end
                j1 = j-1;
                j2 = j;
                if j1 <= 0
                    j1 = length(a2);
                end

               mii(i,:) = [mi(j1),mi(j2),si(i)];
               val1 = (a2(i) - a1(j2)) / (a1(j1) - a1(j2));
               val2 = (a1(j1) - a2(i))/ (a1(j1) - a1(j2));
               valii(i,3) = -1;
               K(si,:)
            end

            for i = 1:length(si)

                j = 1;
                while (a2(i) > a1(j)) && (j <= length(a1))
                    j = j + 1;
                end
                j1 = j-1;
                j2 = j;
                if j1 <= 0
                    j1 = length(a1);
                end

               mii(i,:) = [mi(j1),mi(j2),si(i)];
               val1 = (a2(i) - a1(j2)) / (a1(j1) - a1(j2));
               val2 = (a1(j1) - a2(i))/ (a1(j1) - a1(j2));
               valii(i,3) = -1;
               K(si,:)
            end
            F(si) = 0;

            sii = 1:length(si);
            sii = [sii;sii;sii];
            mii = mii';
            valii = valii';
            K(mi,:) = K(mi,:)
            K(si,:) = sparse(sii(:),mii(:),valii(:),length(si),obj.m.Nn+obj.NcoilArms+obj.Ncoils+obj.NstarConnections);
            

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
            end

            obj.evalBe;
            disp('initial geuss evaluated.')
            toc, disp('-------------------------------------------------------');

            if obj.edata.areAllLinear
                obj.evalHe;
                obj.evalBn;
                obj.evalHn;
                return
            end

            % loop for nonlinear solver
            tic, disp('-------------------------------------------------------');

            % initials values
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
            while (RelError > obj.solverSettings.relativeError) && (Iterations < obj.solverSettings.maxIteration)

                % starting loop time
                loopTime = tic;

                % evaluation of B2 for each elements
                obj.evalBe;
                [obj.solverHistory.totalEnergy(end + 1),obj.solverHistory.totalConergy(end + 1)] = obj.evalTotalEnergyCoenergy;
                Bk = obj.results.Bxg.^2 + obj.results.Byg.^2;

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

                K(si, :) = sparse(1:length(si),si,ones(1,length(si)),length(si),obj.m.Nn+obj.NcoilArms+obj.Ncoils+obj.NstarConnections);
            FF(si) = 0;

                % solving [K][U] = [F]
                dU = K \ FF;
                solVector = full(solVector + alphaNR*dU);

                obj.results.A = solVector(1:obj.m.Nn);
                obj.results.VICoilArms = solVector(obj.m.Nn+1:obj.m.Nn+obj.NcoilArms);
                obj.results.ICoils = solVector(obj.m.Nn+obj.NcoilArms+1:obj.m.Nn+obj.NcoilArms+obj.Ncoils);

                for i = 1:obj.Ncoils
                    % get coil pointer
                    cptr = obj.coils.(coilNames{i});
                    cptr.current(end) = solVector(obj.m.Nn + obj.NcoilArms + cptr.ci);
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
                        alphaNR = min((1.05+2e-2*rand)*alphaNR+2e-3*rand,0.9);
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
            timeHolder = obj.dispLine;

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

            idx = mrptr.idx;
            obj.m.nodes(idx,:) = emdlab_g2d_rotatePoints(obj.m.nodes(idx,:), rotAngle, xc, yc);
            obj.dispMessage('rotation compeleted.', timeHolder);

        end

    end
end
