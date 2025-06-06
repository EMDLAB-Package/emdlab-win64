% developer: https://ComProgExpert.com, Ali Jamali-Fard
% two dimensional nonlinear magnetic-static solver
% nonlinear solver: Newton-Raphson
% first order triangular mesh
% triangular lagrangian elements: 6 points per element
% isotropic
% homogenous
% nonlinear

classdef emdlab_solvers_ms2d_tl6_ihnl < handle & emdlab_solvers_ms2d_tl6
    
    methods
        
        % Initialization
        function obj = emdlab_solvers_ms2d_tl6_ihnl(m)
            
            % mesh pointer
            obj.m = m;
            
            % default settings for solver
            obj.solverSettings.relativeError = 1e-4;
            obj.solverSettings.maxIteration = 20;
            
            % set default properties of mzs
            mzNames = fieldnames(obj.m.mzs);
            
            for i = 1:numel(mzNames)
                obj.setdp(mzNames{i});
            end
            
        end
        
        % Solver
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
        
        function setMonitor(obj, value)
            
            obj.monitorResiduals = value;
            
        end
        
        function assignEdata(obj, InitNur)
            
            % check states
            if obj.isElementDataAssigned, return; end
            
            % preparing mesh data
            obj.m.evalKeMeFe_TL6;
            tic, disp('-------------------------------------------------------');
            
            % assigning material and force data to each triangle
            
            % allocation of memory
            obj.edata.MagneticReluctivity = zeros(3, obj.m.Ne);
            obj.edata.ElectricConductivity = zeros(3, obj.m.Ne);
            obj.edata.InternalCurrentDensity = zeros(3, obj.m.Ne);
            obj.edata.MagnetizationX = zeros(3, obj.m.Ne);
            obj.edata.MagnetizationY = zeros(3, obj.m.Ne);
            
            % getting mesh zones
            mzsName = fieldnames(obj.m.mzs);
            
            % loop over mesh zones
            for i = 1:obj.m.Nmzs
                
                mzptr = obj.m.mzs.(mzsName{i});
                
                if ~obj.m.mts.(mzptr.material).MagneticPermeability.isIsotropic
                    
                    throw(MException('', 'Some materials are non-isotropic, please select the correct solver.'));
                    
                elseif obj.m.mts.(mzptr.material).MagneticPermeability.isLinear
                    
                    % assigning Magnetic Permeability
                    obj.edata.MagneticReluctivity(:,obj.m.ezi(:, mzptr.zi)) = 1/obj.m.mts.(mzptr.material).MagneticPermeability.value;
                    
                else
                    
                    if nargin == 2
                        obj.edata.MagneticReluctivity(:,obj.m.ezi(:, mzptr.zi)) = InitNur * obj.pcts.nu0;
                    else
                        obj.edata.MagneticReluctivity(:,obj.m.ezi(:, mzptr.zi)) = 0.001 * obj.pcts.nu0;
                    end
                    
                end
                
                % assigning Electric Conductivity
                obj.edata.ElectricConductivity(:,obj.m.ezi(:, mzptr.zi)) = obj.m.mts.(mzptr.material).ElectricConductivity.value;
                
                % assigning Internal Current Density
                if mzptr.props.isExcited
                    
                    switch mzptr.props.excitation.type
                        
                        case 'currentDensity'
                            
                            obj.edata.InternalCurrentDensity(:,obj.m.ezi(:, mzptr.zi)) = ...
                                mzptr.props.excitation.value * obj.units.k_currentDensity;
                            
                        case 'current'
                            
                            obj.edata.InternalCurrentDensity(:,obj.m.ezi(:, mzptr.zi)) = ...
                                mzptr.props.excitation.value * obj.units.k_current / (mzptr.getArea * obj.units.k_length^2);
                            
                    end
                    
                end
                
                % assigning Magnetization
                if mzptr.props.isMagnetized
                    
                    M = mzptr.props.magnetization.getM(mzptr.getCenterOfElements);
                    obj.edata.MagnetizationX(:,obj.m.ezi(:, mzptr.zi)) = repmat(M(:, 1)',3,1);
                    obj.edata.MagnetizationY(:,obj.m.ezi(:, mzptr.zi)) = repmat(M(:, 2)',3,1);
                    
                end
                
            end
            
            % applying current of excitation matrices
            mNames = fieldnames(obj.exmtcs);
            
            for i = 1:numel(mNames)
                
                % matrix pointer
                mptr = obj.exmtcs.(mNames{i});

                for j = 1:mptr.Nmzs
                    
                    mzptr = obj.m.mzs.(mptr.mzsName{j});
                    cptr = obj.coils.(mptr.mzsName{j});
                    obj.edata.InternalCurrentDensity(:,obj.m.ezi(:, mzptr.zi)) = ...
                        cptr.sign * cptr.turns * mptr.current * obj.units.k_current  / mptr.np / (mzptr.getArea * obj.units.k_length^2);
                    
                end
                
            end
            
            disp('Initialization of material and force data compeleted.')
            toc, disp('-------------------------------------------------------');
            
            % change states
            obj.isElementDataAssigned = true;
            
        end
        
        function solve(obj)
            
            % prerequisties
            obj.assignEdata;
            
            % updating boundary conditions
            obj.bcs.updateAll;
            
            % getting mesh zone names
            mzNames = fieldnames(obj.m.mzs);
            
            % Construction of [K] and [F]
            tic, disp('-------------------------------------------------------');           

            % assembling the total load vector
            F = (obj.m.mtcs.Fe .* obj.edata.InternalCurrentDensity(1,:)) * obj.units.k_length^2 + ...
                (obj.edata.MagnetizationX(1,:) .* obj.m.mtcs.FeMx + obj.edata.MagnetizationY(1,:) .* obj.m.mtcs.FeMy) * (obj.units.k_length * obj.units.k_magnetisation);
            F = sparse(obj.m.cl', ones(6 * obj.m.Ne, 1), F);
            
            % Assembeling [K]
            [Iindex, Jindex] = getij(6,1);
            Iindex = obj.m.cl(:, Iindex)';
            Jindex = obj.m.cl(:, Jindex)';
            K = sparse(Iindex, Jindex, obj.edata.MagneticReluctivity(1,:).*obj.m.mtcs.Ke1 + obj.edata.MagneticReluctivity(2,:).*obj.m.mtcs.Ke2 + obj.edata.MagneticReluctivity(3,:).*obj.m.mtcs.Ke3);
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
            
            % loop for nonlinear solver
            tic, disp('-------------------------------------------------------');
            
            % initials values
            RelError = inf;
            Iterations = 0;
            xNgt = obj.m.Ne;
            xNgp = obj.m.Nn;
            
            % preparing error monitoring
            if obj.monitorResiduals

                ERF = gcf; cla;
                set(ERF, 'Name', 'md2d_tl6_ihnl solver', 'WindowStyle', 'Normal');
                er = animatedline('color', 'r', 'Linewidth', 1.2, 'Marker', 's', 'MarkerEdgeColor','k');
                title('Relative Error');
                ylabel('log10(||dA||/||A||)');
                xlabel('Iteration Number');
                cAxis = gca;
                box on;
                grid on;
                set(gca, 'ylim', [log10(obj.solverSettings.relativeError)-1,1]);
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

            % updating dnudB2
            dnudB2 = zeros(3, xNgt);

            alphaNR = 0.8;
            
            % loop for non-linearity
            fprintf('Iter|Error   |Residual|time\n');
            while (RelError > obj.solverSettings.relativeError) && (Iterations < obj.solverSettings.maxIteration)
                
                % starting loop time
                loopTime = tic;
                
                % evaluation of B2 for each elements
                obj.evalBe;
                [obj.solverHistory.totalEnergy(end + 1),obj.solverHistory.totalConergy(end + 1)] = obj.evalTotalEnergyCoenergy;
                Bk = sqrt(obj.results.Bxg.^2 + obj.results.Byg.^2);
                
                % updating nu
                for i = 1:obj.m.Nmzs
                    mzptr = obj.m.mzs.(mzNames{i});
                    
                    if ~obj.m.mts.(mzptr.material).MagneticPermeability.isLinear
                        for j = 1:3
                            obj.edata.MagneticReluctivity(j,obj.m.ezi(:, mzptr.zi)) = ...
                                ppval(obj.m.mts.(mzptr.material).vB, ...
                                Bk(j,obj.m.ezi(:, mzptr.zi)));
                        end
                    end
                    
                end
                                
                % updating dnudB2
                for i = 1:obj.m.Nmzs
                    mzptr = obj.m.mzs.(mzNames{i});
                    
                    if ~obj.m.mts.(mzptr.material).MagneticPermeability.isLinear
                        for j = 1:3
                            dnudB2(j,obj.m.ezi(:, mzptr.zi)) = ...
                                ppval(obj.m.mts.(mzptr.material).dvdB, ...
                                Bk(j,obj.m.ezi(:, mzptr.zi))) ./ ...
                                (2 * Bk(j,obj.m.ezi(:, mzptr.zi)));
                        end
                    end
                    
                end
                
                % construction of stiffness matrix [K]
                K = sparse(Iindex, Jindex, obj.edata.MagneticReluctivity(1,:).*obj.m.mtcs.Ke1 + obj.edata.MagneticReluctivity(2,:).*obj.m.mtcs.Ke2 + obj.edata.MagneticReluctivity(3,:).*obj.m.mtcs.Ke3);
                
                % construction of [K] and [F] in NR algorithm
                FF = -K * obj.results.A + F;
                
                % evaluation and adding of jacobian matrix
                [G1, G2, G3] = emdlab_m2d_tl6_evalG(obj.m.cl, obj.m.mtcs.Ke1, obj.m.mtcs.Ke2, obj.m.mtcs.Ke3, obj.m.JIT, obj.results.A, obj.results.Bxg, obj.results.Byg, dnudB2);
                K = K + sparse(Iindex, Jindex,  (G1+G2+G3)/ obj.units.k_length);
                
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

                if alphaNR == 0.8
                    alphaNR = 0.7 + 0.2*2*(rand-0.5);
                else
                    alphaNR = 0.8;
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

            % evaluation of B2 for each elements
            obj.evalBe;
            Bk = sqrt(obj.results.Bxg.^2 + obj.results.Byg.^2);

            % updating nu
            for i = 1:obj.m.Nmzs
                mzptr = obj.m.mzs.(mzNames{i});

                if ~obj.m.mts.(mzptr.material).MagneticPermeability.isLinear
                    for j = 1:3
                        obj.edata.MagneticReluctivity(j,obj.m.ezi(:, mzptr.zi)) = ...
                            ppval(obj.m.mts.(mzptr.material).vB, ...
                            Bk(j,obj.m.ezi(:, mzptr.zi)));
                    end
                end

            end

            obj.isBnEvaluated = false;
            
        end
        
        % Solver History
        function plotRelativeError(obj)
            
            figure('Name', 'emdlab -> ms2d_ihnl solver');
            plot(1:obj.solverHistory.iterations, ...
                log10(obj.solverHistory.relativeError), 'color', 'r', 'Linewidth', 1.2, ...
                'marker', 's', 'MarkerEdgeColor', 'b');
            title('Relative Error (|dA|/|A|)')
            ylabel('log10(|dA|/|A|)')
            xlabel('Iteration Number')
            if obj.solverHistory.iterations == 1, return; end
            set(gca, 'xlim', [1, obj.solverHistory.iterations]);
            
        end
        
        function plotEnergy(obj)

            figure('Name', 'emdlab -> ms2d_ihnl solver');
            plot(1:obj.solverHistory.iterations, ...
                obj.solverHistory.totalEnergy, 'color', 'r', 'Linewidth', 1.2, ...
                'marker', 's', 'MarkerEdgeColor', 'b');
            title('Total Energy')
            ylabel('Energy [J]')
            xlabel('Iteration Number')
            if obj.solverHistory.iterations == 1, return; end
            set(gca, 'xlim', [1, obj.solverHistory.iterations]);

        end
        
        function plotEnergyResidual(obj)

            figure('Name', '[EMDLAB] IHNLNRMSTL3 Solver', 'NumberTitle', 'off');
            plot(2:obj.solverHistory.iterations, ...
                diff(obj.solverHistory.totalEnergy), 'color', 'r', 'Linewidth', 1.2, ...
                'marker', 's', 'MarkerEdgeColor', 'b');
            title('Residual of Energy')
            ylabel('Residual [J]')
            xlabel('Iteration Number')
            if obj.solverHistory.iterations == 1, return; end
            set(gca, 'xlim', [1, obj.solverHistory.iterations]);
            
        end
        
    end
    
end
