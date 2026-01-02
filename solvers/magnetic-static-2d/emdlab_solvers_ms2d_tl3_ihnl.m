% EMDLAB: Electrical Machines Design Laboratory
% two dimensional nonlinear magnetic-static solver
% nonlinear solver: Newton-Raphson
% first order triangular mesh
% triangular lagrangian elements: 3 points per element
% isotropic
% homogenous
% nonlinear

classdef emdlab_solvers_ms2d_tl3_ihnl < handle & emdlab_solvers_ms2d_tlcp

    methods

        % initialization
        function obj = emdlab_solvers_ms2d_tl3_ihnl(m)

            % mesh pointer
            m.ggmesh;
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
            obj.m.evalKeMeFe_TL3;
            tic, disp('-------------------------------------------------------');

            % assigning material and force data to each triangle

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

                % get pointer to mesh zone
                mzptr = obj.m.mzs.(mzsName{i});

                % assigning magnetic reluctivity
                if ~obj.m.mts.(mzptr.material).MagneticPermeability.isIsotropic

                    throw(MException('', 'Some materials are non-isotropic, please select the correct solver.'));

                elseif obj.m.mts.(mzptr.material).MagneticPermeability.isLinear

                    obj.edata.MagneticReluctivity(obj.m.ezi(:, mzptr.zi)) = 1/obj.m.mts.(mzptr.material).MagneticPermeability.value;

                else

                    obj.edata.areAllLinear = false;
                    if nargin == 2
                        obj.edata.MagneticReluctivity(obj.m.ezi(:, mzptr.zi)) = obj.pcts.nu0 / initialRelativePermeability;
                    else
                        obj.edata.MagneticReluctivity(obj.m.ezi(:, mzptr.zi)) = obj.pcts.nu0 / 500;
                    end

                end

                % assigning electric conductivity
                obj.edata.ElectricConductivity(obj.m.ezi(:, mzptr.zi)) = obj.m.mts.(mzptr.material).ElectricConductivity.value;

                % assigning internal current density
                if mzptr.props.isExcited

                    switch mzptr.props.excitation.type

                        case 'currentDensity'

                            obj.edata.InternalCurrentDensity(obj.m.ezi(:, mzptr.zi)) = ...
                                mzptr.props.excitation.value * obj.units.k_currentDensity;

                        case 'current'

                            obj.edata.InternalCurrentDensity(obj.m.ezi(:, mzptr.zi)) = ...
                                mzptr.props.excitation.value * obj.units.k_current / (mzptr.getArea * obj.units.k_length^2);

                    end

                end

                % assigning magnetization
                if mzptr.props.isMagnetized

                    M = mzptr.props.magnetization.getM(mzptr.getCenterOfElements);
                    obj.edata.MagnetizationX(obj.m.ezi(:, mzptr.zi)) = M(:, 1)';
                    obj.edata.MagnetizationY(obj.m.ezi(:, mzptr.zi)) = M(:, 2)';

                end

            end

            % applying current of excitation matrices
            coilNames = fieldnames(obj.coils);

            for i = 1:numel(coilNames)

                % get coil pointer
                cptr = obj.coils.(coilNames{i});

                for j = 1:cptr.NcoilArms

                    % get coil arm pointer
                    mzptr = obj.m.mzs.(cptr.coilArms(j));

                    % set coil arm current density
                    obj.edata.InternalCurrentDensity(obj.m.ezi(:, mzptr.zi)) = ...
                        mzptr.props.direction * mzptr.props.turns * cptr.current * obj.units.k_current ...
                        / (mzptr.getArea * obj.units.k_length^2);

                end

            end

            disp('Initialization of material and force data compeleted.')
            toc, disp('-------------------------------------------------------');

            % change states
            obj.isElementDataAssigned = true;

        end

        % solver core
        function solve(obj, varargin)

            % prerequisties
            obj.assignEdata(varargin{:});

            % updating boundary conditions
            obj.bcs.updateAll;

            % Construction of [K] and [F]
            tic, disp('-------------------------------------------------------');

            % Assembeling [F]
            % assembling the total load vector
            F = (obj.m.mtcs.Fe .* obj.edata.InternalCurrentDensity) * obj.units.k_length^2 + ...
                (obj.edata.MagnetizationX .* obj.m.mtcs.FeMx + obj.edata.MagnetizationY .* obj.m.mtcs.FeMy) * (obj.units.k_length * obj.units.k_magnetisation);
            F = sparse(obj.m.cl', ones(3 * obj.m.Ne, 1), F);

            % Assembeling [K]
            [Iindex, Jindex] = emdlab_flib_getij(3,1);
            Iindex = obj.m.cl(:, Iindex)';
            Jindex = obj.m.cl(:, Jindex)';
            K = sparse(Iindex, Jindex, obj.edata.MagneticReluctivity .* obj.m.mtcs.Ke);
            disp('Construction of [K] and [F] compeleted.');
            toc, disp('-------------------------------------------------------');
            tic, disp('-------------------------------------------------------');

            % imposing boundary conditions on [K] and [F]
            % dbcs
            if obj.bcs.Nd
                F(obj.bcs.iD) = obj.bcs.vD;
                K(obj.bcs.iD, :) = sparse(1:obj.bcs.Ndbcs, obj.bcs.iD, ones(1, obj.bcs.Ndbcs), obj.bcs.Ndbcs, obj.m.Nn);
            end

            % opbcs
            if obj.bcs.Nop
                F(obj.bcs.mOP) = F(obj.bcs.mOP) - F(obj.bcs.sOP);
                F(obj.bcs.sOP) = 0;
                K(obj.bcs.mOP, :) = K(obj.bcs.mOP, :) - K(obj.bcs.sOP, :);
                K(obj.bcs.sOP, :) = sparse([1:obj.bcs.Nopbcs, 1:obj.bcs.Nopbcs], ...
                    [obj.bcs.mOP; obj.bcs.sOP], ones(1, 2 * obj.bcs.Nopbcs), obj.bcs.Nopbcs, obj.m.Nn);
            end

            % epbcs
            if obj.bcs.Nep
                F(obj.bcs.mEP) = F(obj.bcs.mEP) + F(obj.bcs.sEP);
                F(obj.bcs.sEP) = 0;
                K(obj.bcs.mEP, :) = K(obj.bcs.mEP, :) + K(obj.bcs.sEP, :);
                K(obj.bcs.sEP, :) = sparse([1:obj.bcs.Nepbcs, 1:obj.bcs.Nepbcs], ...
                    [obj.bcs.mEP; obj.bcs.sEP], [ones(1, obj.bcs.Nepbcs), -ones(1, obj.bcs.Nepbcs)], obj.bcs.Nepbcs, obj.m.Nn);
            end

            disp('All boundary condition imposed.');
            toc, disp('-------------------------------------------------------');

            % solving [K][U] = [F]
            tic, disp('-------------------------------------------------------');

            % solving equation KU = F
            if ~any(F)
                obj.results.A = full(F);
                return
            end

            obj.results.A = full(K \ F);
            obj.evalBe;
            disp('initial geuss evaluated.')
            toc, disp('-------------------------------------------------------');
            if obj.edata.areAllLinear
                obj.evalHe;
                obj.evalBn;
                obj.evalHn;
                [obj.solverHistory.totalEnergy,obj.solverHistory.totalConergy] = obj.evalTotalEnergyCoenergy;
                obj.solverHistory.relativeError = 0;
                obj.solverHistory.iterations = 1;
                obj.solverHistory.energyResidual = [];
                return
            end

            % loop for nonlinear solver
            tic, disp('-------------------------------------------------------');

            % initials values
            RelEResidual = inf;
            RelError = inf;
            Iterations = 0;
            xNgt = obj.m.Ne;
            xNgp = obj.m.Nn;

            % preparing error monitoring
            if obj.monitorResiduals

                ERF = gcf; cla;
                set(ERF, 'Name', 'md2d_tl3_ihnl solver', 'WindowStyle', 'Normal');
                er = animatedline('color', 'r', 'Linewidth', 1.2, 'Marker', 's', 'MarkerEdgeColor','k');
                title('Relative Error');
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
                K = sparse(Iindex, Jindex, obj.edata.MagneticReluctivity .* obj.m.mtcs.Ke);

                % construction of [K] and [F] in NR algorithm
                FF = -K * obj.results.A + F;

                % evaluation and adding of jacobian matrix
                K = K + sparse(Iindex, Jindex, emdlab_m2d_tl3_evalG(obj.m.cl, obj.m.mtcs.Ke, obj.m.JIT, obj.results.A, dnudB2) / obj.units.k_length^2);

                % imposing boundary conditions on incrimentals
                % dbcs
                if obj.bcs.Nd
                    FF(obj.bcs.iD) = 0;
                    K(obj.bcs.iD, :) = sparse(1:obj.bcs.Ndbcs, obj.bcs.iD, 1, obj.bcs.Ndbcs, xNgp);
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
                obj.results.A = full(obj.results.A + alphaNR*dU);

                % check for convergency
                Residual = norm(dU, 2);
                RelError = Residual / norm(obj.results.A, 2);

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

            % update magnetic flux density
            obj.evalBe;
            obj.evalBn;

            % update magnetic reluctivity
            Bk = obj.results.Bxg.^2 + obj.results.Byg.^2;

            % updating nu -> B2
            for i = 1:obj.m.Nmts
                mtptr = obj.m.mts.(obj.m.materialNames(i));

                if ~mtptr.MagneticPermeability.isLinear
                    obj.edata.MagneticReluctivity(obj.m.emi(i,:)) = ppval(mtptr.vB2, Bk(obj.m.emi(i,:)));
                end

            end

            % update magnetic flux intensity
            obj.evalHe;
            obj.evalHn;

            % change state
            obj.isResultsValid = true;

        end

        % solver history plots
        function varargout = plotRelativeError(obj, varargin)

            [f,ax] = emdlab_flib_fax(varargin{:});
            title(ax, 'emdlab -> ms2d_ihnl solver');
            semilogy(1:obj.solverHistory.iterations, ...
                obj.solverHistory.relativeError, 'color', 'r', 'Linewidth', 1.2, ...
                'marker', 's', 'MarkerEdgeColor', 'b', 'parent', ax);
            title('Relative Error (|dA|/|A|)')
            ylabel('|dA|/|A|')
            xlabel('Iteration Number')
            if obj.solverHistory.iterations == 1, return; end
            set(gca, 'xlim', [1, obj.solverHistory.iterations]);
            grid(ax,'on');

            if nargout == 1, varargout{1} = f;
            elseif nargout == 2, varargout{1} = f; varargout{2} = ax;
            elseif nargout > 1, error('Too many output argument.');
            end

        end

        function varargout = plotEnergy(obj, varargin)

            [f,ax] = emdlab_flib_fax(varargin{:});
            title('emdlab -> ms2d_ihnl solver');
            plot(1:obj.solverHistory.iterations, ...
                obj.solverHistory.totalEnergy, 'color', 'r', 'Linewidth', 1.2, ...
                'marker', 's', 'MarkerEdgeColor', 'b', 'parent', ax);
            ylabel('Energy [J]')
            xlabel('Iteration Number')
            legend('Total Energy = ' + string(obj.solverHistory.totalEnergy(end)) + 'J', 'FontSize', 14)
            if obj.solverHistory.iterations == 1, return; end
            set(ax, 'xlim', [1, obj.solverHistory.iterations]);
            grid(ax,'on');

            if nargout == 1, varargout{1} = f;
            elseif nargout == 2, varargout{1} = f; varargout{2} = ax;
            elseif nargout > 1, error('Too many output argument.');
            end

        end

        function varargout = plotCoenergy(obj, varargin)

            [f,ax] = emdlab_flib_fax(varargin{:});
            title(ax,'emdlab -> ms2d_ihnl solver');
            plot(1:obj.solverHistory.iterations, ...
                obj.solverHistory.totalConergy, 'color', 'r', 'Linewidth', 1.2, ...
                'marker', 's', 'MarkerEdgeColor', 'b', 'parent', ax);
            ylabel('Coenergy [J]')
            xlabel('Iteration Number')
            legend('Total Coenergy = ' + string(obj.solverHistory.totalConergy(end)) + 'J', 'FontSize', 14)
            if obj.solverHistory.iterations == 1, return; end
            set(ax, 'xlim', [1, obj.solverHistory.iterations]);
            grid(ax,'on');

            if nargout == 1, varargout{1} = f;
            elseif nargout == 2, varargout{1} = f; varargout{2} = ax;
            elseif nargout > 1, error('Too many output argument.');
            end

        end

        function varargout = plotEnergyResidual(obj, varargin)

            [f,ax] = emdlab_flib_fax(varargin{:});
            title(ax,'emdlab -> ms2d_ihnl solver');
            relERERR = abs(diff(obj.solverHistory.totalEnergy))./obj.solverHistory.totalEnergy(2:end);
            semilogy(2:obj.solverHistory.iterations, ...
                relERERR, 'color', 'r', 'Linewidth', 1.2, ...
                'marker', 's', 'MarkerEdgeColor', 'b');
            title('Relative Energy Residual')
            ylabel('Relative Energy Residual')
            xlabel('Iteration Number')
            if obj.solverHistory.iterations == 1, return; end
            set(gca, 'xlim', [1, obj.solverHistory.iterations]);
            grid(ax,'on');

            if nargout == 1, varargout{1} = f;
            elseif nargout == 2, varargout{1} = f; varargout{2} = ax;
            elseif nargout > 1, error('Too many output argument.');
            end

        end

        function gui(obj)

            f = gui@emdlab_solvers_ms2d_tlcp(obj);

            tpg = findobj(f, 'Tag', 'MainTabGroup');

            t = uitab(tpg,'Title','Solver Performance');
            tgp_sub = uitabgroup(t);

            t = uitab(tgp_sub,'Title','Energy');
            ax = axes(t);
            obj.plotEnergy(ax);

            t = uitab(tgp_sub,'Title','Conergy');
            ax = axes(t);
            obj.plotCoenergy(ax);

            t = uitab(tgp_sub,'Title','Relative Error');
            ax = axes(t);
            obj.plotRelativeError(ax);

            t = uitab(tgp_sub,'Title','Relative Energy Residual');
            ax = axes(t);
            obj.plotEnergyResidual(ax);

            t = uitab(tpg,'Title','About');
            uicontrol('Parent', t, ...
                'Style', 'edit', ...
                'Units', 'normalized', ...
                'Position', [0.05 0.05 0.9 0.9], ...
                'String', { ...
                'EMDLAB', ...
                'Electrical Machines Design Laboratory', '', ...
                'Repository link:', ...
                'https://github.com/EMDLAB-Package/emdlab-win64', '', ...
                'Paper link:', ...
                'https://www.sciencedirect.com/science/article/pii/S2352711025004121'} , ...
                'HorizontalAlignment', 'center', ...
                'FontSize', 14, ...
                'Max', 2);

            set(f,'visible', 'on');

        end

    end

end
