classdef emdlab_solvers_opt_nsga2 < handle

    properties

        nPop = 5; % population size in each iteration
        indexPop (1,:) double {mustBeInteger, mustBePositive}; % index of population members in genoTypes

        pc = 40; % crossover percentage
        pm = 30; % mutation percentage

        nFails = 0; % number of fails for function evaluations
        nfe = 0; % number of function evaluations
        nfe_max = 1e4; % maximum number of function evaluations
        acceptableCost = 1e-2; % acceptable cost to stop the solver

        gamma = 10; % extrapolation parameter used in the arichmetic crossover
        Nmu = 1; % number of gens used for mutation

        % matrices to store full information
        genoTypes (:,:) double;
        fenoTypes (:,:) double;
        costs (:,:) double;
        bestCosts (:,:) double;
        bestIndividuals (:,1) double;       

        costFunction; % cost function handle

        conFlag = true; % if true solver works in the specefied domain only

        % how to treat constraints
        deathPenaltyFlag = true;
        penaltyFunctionFlag = false;

    end

    properties (Hidden=true)

        Nx; % number of genoTypes
        Ny; % number of fenoTypes
        Nobj; % number of objectives: 1 for single objective, >1 for multi-objective
        VarMin; % minimum value of variables
        VarMax; % maximum value of variables

    end

    properties (Dependent=true)

        nc; % number of parents used for crossover = number of offsprings
        nm; % number of individuals used for mutation

    end

    methods
        %% constructor
        function obj = emdlab_solvers_opt_nsga2(costFunction, VarMin, VarMax)

            % check function handle
            if ~isa(costFunction, 'function_handle')
                error('Please input a function handle.')
            end
            obj.costFunction = costFunction;

            % check varMin and varMax
            if length(VarMin) ~= length(VarMax)
                error('The size of varMin vector must be equal to the size of varMax vector.')
            end

            if any(VarMax < VarMin)
                error('All elements of varMax vector must be higher than their corresponding elements in varMin vector.');
            end

            % run the function handle once to check if it is working properly
            try
                [costVector, outputs, ~] = obj.costFunction((VarMin+VarMax)/2);                
            catch
                error('Cost function is not working properly. Make sure that design is feasible and cost function returns [costs, fenoTypes, gx].');
            end

            % set Nx & Ny & Nobj
            obj.VarMin = VarMin;
            obj.VarMax = VarMax;
            obj.Nx = numel(VarMin);
            obj.Ny = numel(outputs);
            obj.Nobj = length(costVector);

            % initialize inputs/outputs
            obj.genoTypes = zeros([],obj.Nx);
            obj.fenoTypes = zeros([],obj.Ny);
            obj.costs = zeros([],obj.Nobj);

        end

        %% get methods
        function y = get.nc(obj)
            y = 2*floor(obj.pc*1e-2*obj.nPop/4);
            y = min(y,floor(obj.nPop/4));
        end

        function y = get.nm(obj)
            y = round(obj.pm*1e-2*obj.nPop);
        end

        %% set methods
        function setPopulationSize(obj, value)
            obj.nPop = value;
        end

        function setCrossoverPercentage(obj, value)
            obj.pc = value;
        end

        function setCrossoverExtrapolationPercentage(obj, value)
            obj.gamma = value;
        end

        function setMutationPercentage(obj, value)
            obj.pm = value;
        end

        function setAcceptableCost(obj, value)
            obj.acceptableCost = value;
        end

        function setMutationRate(obj, value)
            obj.Nmu = ceil(value * 1e-2 * obj.Nx);
        end

        function enableBoundedSreach(obj)
            obj.conFlag = true;
        end

        function disableBoundedSreach(obj)
            obj.conFlag = false;
        end

        function setMaxFunctionEvaluations(obj,value)
            obj.nfe_max = value;
        end

        %% initialization
        function generateInitialPopulationAndEvaluate(obj)

            % generate initial genoTypes randomly
            obj.genoTypes = obj.VarMin + rand(obj.nPop,obj.Nx).*(obj.VarMax - obj.VarMin);

            % allocate memory from RAM to extend storage matrices
            obj.fenoTypes = zeros(obj.nPop,obj.Ny);
            obj.costs = zeros(obj.nPop,obj.Nobj);

            % index vector pointing to live population
            obj.indexPop = 1:obj.nPop;

            for i = obj.indexPop
                % evaluate cost function
                try
                    [obj.costs(i,:), obj.fenoTypes(i,:), ~] = obj.costFunction(obj.genoTypes(i,:));
                catch
                    error('Unfeasible designs in initial population. Check variable bound for feasible designs.');
                end
            end
            obj.nfe = obj.nfe + obj.nPop;

            % sort vs. cost
            if obj.Nobj == 1
                [obj.costs, sortOrder] = sort(obj.costs);
            else
                [obj.costs, sortOrder] = nd_sort(obj.costs);
            end
            obj.genoTypes(sortOrder,:) = obj.genoTypes(sortOrder,:);
            obj.fenoTypes(sortOrder,:) = obj.fenoTypes(sortOrder,:);

            % store the best cost
            obj.bestCosts(end+1,:) = obj.costs(1,:);
            obj.bestIndividuals(end+1) = 1;

        end

        function generateInitialPopulationAndEvaluateInParallel(obj)

            % generate initial genoTypes randomly
            obj.genoTypes = obj.VarMin + rand(obj.nPop,obj.Nx).*(obj.VarMax - obj.VarMin);

            % index vector pointing to live population
            obj.indexPop = 1:obj.nPop;

            % define proper matrices for parfor loop
            genoTypesInParfor = obj.genoTypes;
            fenoTypesInParfor = zeros(obj.nPop,obj.Ny);
            costsInParfor = zeros(obj.nPop,obj.Nobj);
            localCostFunction = @(x) obj.costFunction(x);

            % parfor loop to evalue initial population
            parfor i = obj.indexPop
                % evaluate cost function
                try 
                    [costsInParfor(i,:), fenoTypesInParfor(i,:), ~] = localCostFunction(genoTypesInParfor(i,:));
                catch
                    error('Unfeasible designs in initial population.');
                end
            end
            obj.nfe = obj.nfe + obj.nPop;

            % allocate memory from RAM to extend storage matrices
            obj.fenoTypes = fenoTypesInParfor;
            obj.costs = costsInParfor;

            % sort vs. cost
            if obj.Nobj == 1
                [obj.costs, sortOrder] = sort(obj.costs);
            else
                [obj.costs, sortOrder] = nd_sort(obj.costs);
            end
            obj.genoTypes(sortOrder,:) = obj.genoTypes(sortOrder,:);
            obj.fenoTypes(sortOrder,:) = obj.fenoTypes(sortOrder,:);

            % store the best cost
            obj.bestCosts(end+1) = obj.costs(1);
            obj.bestIndividuals(end+1) = 1;

        end

        function generateInitialFeasiblePopulationAndEvaluate(obj, suggestion)

            % allocate memory from RAM to extend storage matrices
            obj.genoTypes = zeros(obj.nPop,obj.Nx);
            obj.fenoTypes = zeros(obj.nPop,obj.Ny);
            obj.costs = zeros(obj.nPop,obj.Nobj);

            % index vector pointing to live population
            obj.indexPop = 1:obj.nPop;

            % loop to generate initial & feasible population            
            for i = obj.indexPop
                while true

                    if i == 1
                        if nargin == 2
                            obj.genoTypes(i,:) = suggestion;
                        else
                            % generate an initial genoTypes vector randomly
                            obj.genoTypes(i,:) = obj.VarMin + rand(1,obj.Nx).*(obj.VarMax - obj.VarMin);
                        end
                    else
                        obj.genoTypes(i,:) = obj.mutateIndividuals(i-1);
                    end

                    % evaluate cost function
                    try
                        [obj.costs(i,:),obj.fenoTypes(i,:), gx] = obj.costFunction(obj.genoTypes(i,:));
                        obj.nfe = obj.nfe + 1;
                        if all(gx<=0), break; end
                        if (i == 1) && (nargin == 2)
                            warning('Your suggestion is not feasible. The solver is trying to find a new feasible one.');  
                            uDir = rand(1,obj.Nx);
                            uDir = uDir / norm(uDir);
                            uAMP = 0.1*(obj.VarMax - obj.VarMin);
                            suggestion = suggestion + uAMP.*uDir;
                        end
                    catch
                        obj.nFails = obj.nFails + 1;
                        if obj.nFails > 50
                            error('Too many fails. The problem is ill conditioned.');
                        end
                    end
                    
                end
            end

            % sort vs. cost
            if obj.Nobj == 1
                [obj.costs, sortOrder] = sort(obj.costs);
            else
                [obj.costs, sortOrder] = nd_sort(obj.costs);
            end
            obj.genoTypes(sortOrder,:) = obj.genoTypes(sortOrder,:);
            obj.fenoTypes(sortOrder,:) = obj.fenoTypes(sortOrder,:);

            % store the best cost
            obj.bestCosts(end+1) = obj.costs(1);
            obj.bestIndividuals(end+1) = 1;

        end

        %% selections: parents for cross over or individuals for mutation
        % parents might be non-unique
        function parents = selectParentsRandomly(obj)

            % select pairs randomly
            index1 = randi([1,obj.nPop],obj.nc,1);
            index2 = randi([1,obj.nPop],obj.nc,1);
            parents = [obj.indexPop(index1);obj.indexPop(index2)];

        end

        % select unique parents 
        function parents = selectUniqueParentsRandomly(obj)

            % allocate memory from RAM
            index1 = zeros(1,obj.nc);
            index2 = zeros(1,obj.nc);

            index = 1:obj.nPop;
            for i = 1:obj.nc

                index1(i) = index(randi(length(index)));
                index = setdiff(index,index1(i));

                index2(i) = index(randi(length(index)));
                index = setdiff(index,index2(i));

            end

            parents = [obj.indexPop(index1);obj.indexPop(index2)];

        end

        % tournament selection
        function parents = selectParentsByRouletteWheel(obj, beta)

            % initialize parent matrix
            parents = zeros(2,obj.nc);

            % calculate probabilities
            p = exp(-beta*obj.costs(obj.indexPop)/max(obj.costs(obj.indexPop)));

            % normalize probabilities
            p = p/sum(p);
            cp = cumsum(p);

            % select pairs by roulette wheel
            for i = 1:obj.nc
                index1 = find(rand<=cp,1,'first');
                index2 = find(rand<=cp,1,'first');
                parents(:,i) = [obj.indexPop(index1);obj.indexPop(index2)];
            end

        end

        % roulette wheel selection
        function parents = selectParentsByTournament(obj, m)

            % initialize parent matrix
            parents = zeros(2,obj.nc);

            % select pairs by tournament
            for i = 1:obj.nc
                index1 = obj.indexPop(randperm(obj.nPop, m));
                index2 = obj.indexPop(randperm(obj.nPop, m));

                if obj.Nobj == 1
                    % tournament for single-opjective
                    [~,index1] = min(obj.costs(index1));
                    [~,index2] = min(obj.costs(index2));
                else
                    % tournament for multi-opjective
                    [~,index1] = nd_sort(obj.costs(index1,:));
                    index1 = index1(1);
                    [~,index2] = nd_sort(obj.costs(index2,:));
                    index2 = index2(1);
                end

                parents(:,i) = [obj.indexPop(index1);obj.indexPop(index2)];
            end

        end

        % individuals might be non-unique
        function individuals = selectIndividualsRandomly(obj)

            % random individuals
            individuals = obj.indexPop(randi([1,obj.nPop],obj.nm,1));

        end

        function individuals = selectUniqueIndividualsRandomly(obj)

            % allocate memory from RAM
            individuals = zeros(1,obj.nm);

            index = 1:obj.nPop;
            for i = 1:obj.nm

                individuals(i) = index(randi(length(index)));
                index = setdiff(index,individuals(i));

            end

            individuals = obj.indexPop(individuals);

        end

        %% new generation
        function offsprings = crossoverParents(obj, parents)

            % allocate memory from RAM
            offsprings = zeros(obj.nc, obj.Nx);

            GAMMA = obj.gamma*1e-2;
            for i = 1:2:obj.nc

                alpha = -GAMMA + (1+2*GAMMA) * rand(1,obj.Nx);

                % child #1
                offsprings(i,:) = alpha.*obj.genoTypes(parents(1,i),:)+(1-alpha).*obj.genoTypes(parents(2,i),:);

                % child #2
                offsprings(i+1,:) = alpha.*obj.genoTypes(parents(2,i),:)+(1-alpha).*obj.genoTypes(parents(1,i),:);

            end

            if obj.conFlag
                offsprings = max(offsprings,obj.VarMin);
                offsprings = min(offsprings,obj.VarMax);
            end

        end

        function mutations = mutateIndividuals(obj, individuals)

            % allocate memory from RAM
            mutations = zeros(length(individuals), obj.Nx);

            sigma = 0.1 * (obj.VarMax-obj.VarMin);
            for i = 1:length(individuals)

                % copy chromosome
                mutations(i,:) = obj.genoTypes(individuals(i),:);

                % select gens randomly to mutate
                j = randperm(obj.Nx,obj.Nmu);

                % mutate selected gens
                mutations(i,j) = mutations(i,j) + sigma(j) .* randn(1,obj.Nmu);

            end

            if obj.conFlag
                mutations = max(mutations,obj.VarMin);
                mutations = min(mutations,obj.VarMax);
            end

        end

        %% evaluate new generation
        function evaluetMergeSortTruncate(obj, offsprings, mutations)
           
            % size of new generation
            NnewPop = size(offsprings,1) + size(mutations,1); 

            % indicies refering to individuals of new generation
            indexNewPop = (1:NnewPop) + size(obj.genoTypes,1);

            % allocate memory from RAM and extend storage matrices
            obj.genoTypes = [obj.genoTypes;[offsprings;mutations]];
            obj.fenoTypes = [obj.fenoTypes;zeros(NnewPop,obj.Ny)];
            obj.costs = [obj.costs;zeros(NnewPop,obj.Nobj)];

            % for loop for evaluation
            for i = indexNewPop
                % evaluate cost function
                try
                    [obj.costs(i,:),obj.fenoTypes(i,:),gx] = obj.costFunction(obj.genoTypes(i,:));
                    obj.nfe = obj.nfe + 1;
                    if any(gx>0), obj.costs(i,:) = inf; end
                catch
                    any(obj.costs(i,:) == inf);
                    obj.nFails = obj.nFails + 1;
                    if obj.nFails > 50
                        error('Too many fails. The problem is ill conditioned.');
                    end
                end
            end            

            % merge: index of live population
            indexLivePop = [obj.indexPop, indexNewPop];

            % sort evaluated live population
            if obj.Nobj == 1
                [~, indexSort] = sort(obj.costs(indexLivePop,:));
            else
                [~, indexSort] = nd_sort(obj.costs(indexLivePop,:));
            end            
            indexLivePop = indexLivePop(indexSort);

            % save best cost
            obj.bestCosts(end+1,:) = obj.costs(indexLivePop(1),:);

            % truncate
            obj.indexPop = indexLivePop(1:obj.nPop);
            obj.bestIndividuals(end+1) = indexLivePop(1);

        end

        function evaluetInparallelMergeSortTruncate(obj, offsprings, mutations)

            % evalue new population
            NnewPop = size(offsprings,1) + size(mutations,1); % size of new generation
            indexNewPop = (1:NnewPop) + size(obj.genoTypes,1);

            % allocate memory from RAM and extend storage matrices
            obj.genoTypes = [obj.genoTypes;[offsprings;mutations]];

            % define propermatrices for parfor loop
            genoTypesInParfor = [offsprings;mutations];
            fenoTypesInParfor = zeros(NnewPop,obj.Ny);
            costsInParfor = zeros(NnewPop,obj.Nobj);
            nFailsInParfor = zeros(1,NnewPop);
            localCostFunction = @(x) obj.costFunction(x);

            % parfor loop for evaluation
            parfor i = 1:NnewPop
                % evaluate cost function
                try
                    [costsInParfor(i,:),fenoTypesInParfor(i,:),gx] = localCostFunction(genoTypesInParfor(i,:));                    
                    if any(gx>0), costsInParfor(i,:) = inf; end
                catch
                    costsInParfor(i,:) = inf;
                    nFailsInParfor(i) = 1;                    
                end
            end    
            obj.nfe = obj.nfe + NnewPop;
            obj.nFails = sum(nFailsInParfor);
            if obj.nFails > 50
                error('Too many fails. The problem is ill conditioned.');
            end

            % allocate memory from RAM and extend storage matrices
            obj.fenoTypes = [obj.fenoTypes;fenoTypesInParfor];
            obj.costs = [obj.costs;costsInParfor];

            % index of live population
            indexLivePop = [obj.indexPop, indexNewPop];

            % sort evaluated live population
            if obj.Nobj == 1
                [~, indexSort] = sort(obj.costs(indexLivePop,:));
            else
                [~, indexSort] = nd_sort(obj.costs(indexLivePop,:));
            end            
            indexLivePop = indexLivePop(indexSort);

            % save best cost
            obj.bestCosts(end+1) = obj.costs(indexLivePop(1));

            % truncate
            obj.indexPop = indexLivePop(1:obj.nPop);
            obj.bestIndividuals(end+1) = indexLivePop(1);

        end

        function trFlag = checkForTermination(obj)

            trFlag = false;
            if obj.Nobj == 1
                if (obj.bestCosts(end) <= obj.acceptableCost) || (obj.nfe >= obj.nfe_max)
                    trFlag = true;
                end
            else
                if (obj.nfe >= obj.nfe_max)
                    trFlag = true;
                end
            end

        end

        %% visualization functions
        function plotBestCostsVSIteration(obj, logFlag)

            if nargin<2, logFlag = false; end

            figure; box on;
            if logFlag
                semilogy(1:length(obj.bestCosts), obj.bestCosts, 'LineWidth', 2);
            else
                plot(1:length(obj.bestCosts), obj.bestCosts, 'LineWidth', 2);
            end
            grid on;
            xlabel('Iteration');
            ylabel('Cost');
            title("NFE = " + string(obj.nfe) + " Best Cost = " + sprintf('%.2g',obj.bestCosts(end)));

        end

        function printBestIndividual(obj)
            disp('costs:');
            disp(obj.costs(obj.bestIndividuals(end),:));
            disp('genoType:');
            disp(obj.genoTypes(obj.bestIndividuals(end),:));
            disp('fenoTypes:');
            disp(obj.fenoTypes(obj.bestIndividuals(end),:));
        end

        function exportLastPopulation(obj)
            writematrix([obj.genoTypes(obj.indexPop,:),obj.fenoTypes(obj.indexPop,:)],'lastPopulation.csv');
        end

    end

    methods (Static=true)

        function [sortOrder, fronts, ranks, crowd] = ndsort(costsVector)



        end
    end

end