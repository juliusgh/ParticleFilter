classdef ParticleFilter < handle
    %PARTICLEFILTER Particle Filter for Set-Membership Filtering
    %   Contains various algorithms and the related data structures
    
    properties
        system % dynamic system
        algorithm % PFAlgorithms enumeration
        options % options for the Particle Filter
        particles % array that contains all particles
        counts % particle counts
        time % execution time
        areas
        sets
        t
        maxk
        intvS
        intvE
        intvN
        intvSz
        intvEz
        intvNz
    end
    
    methods
        function obj = ParticleFilter(system,algorithm,options)
            %PARTICLEFILTER Construct an instance of this class
            %   Constructs Particle Filter for a system with given options
            if nargin < 3
                options = PFOptions();
                if nargin < 2
                    algorithm = PFAlgorithms.Naive;
                end
            end
            obj.system = system;
            obj.algorithm = algorithm;
            obj.options = options;
            obj.initialize();
        end
        
        function obj = initialize(obj)
            %INITIALIZE Summary of this method goes here
            %   Detailed explanation goes here
            if obj.options.simulate
                obj.system = obj.system.simulate(obj.options.iterations);
            else
%                  obj.system.x = obj.options.xvals;
%                  obj.system.z = obj.options.zvals;
            end
            obj.particles = NaN(obj.system.dim,obj.options.samples,0);
            obj.counts = zeros(1,obj.options.iterations);
            obj.time = NaN;
            obj.t = 0;
            obj.maxk = 0;
        end
        
        function obj = run(obj, divisions)
            %RUN Execute the Particle Filter
            %   Execution time is stored in property time
            if nargin < 2
                divisions = 1;
            end
            obj.initialize();
            if divisions > 1
                dtime = tic;
                subsystems = obj.system.divide(divisions);
                n = numel(subsystems);
                obj.areas = NaN(divisions);
                subPF = ParticleFilter.empty(0,n);
                suboptions = obj.options;
                suboptions.samples = ceil(obj.options.samples / n);
                for i = 1:n % divide
                    subPF(i) = ParticleFilter(subsystems(i),obj.algorithm,suboptions);
                end
                dividetime = toc(dtime);
                if obj.options.debug
                    fprintf('--- finished divide after %f seconds\n',dividetime);
                end
                ctime = tic;
                p = 0;
                par = obj.hasParallel();
                if par
                    D = parallel.pool.DataQueue;
                    if obj.options.debug
                        afterEach(D, @updateStatus);
                    end
                end
                iter = obj.options.iterations;
                parfor i = 1:n % conquer
                    subPF(i) = subPF(i).run();
                    a = subPF(i).area(iter);
                    if par
                        send(D,[i,a]);
                    end
                end
                conquertime = toc(ctime);
                if obj.options.debug
                    fprintf('--- finished conquer after %f seconds\n',conquertime);
                end
                mtime = tic;
%                 for i = 1:n
%                     obj.areas(i) = subPF(i).area(obj.options.iterations);
%                 end
                % start with empty particle array
                obj.particles = NaN(obj.system.dim,0,0);
                % merge smallest particle arrays first
                [~,ind] = sort([subPF.maxk]);
                subPF = subPF(ind);
                for i = 1:n % merge
                    obj.merge(subPF(i));
                    subPF(i).clear();
                end
                mergetime = toc(mtime);
                if obj.options.debug
                    fprintf('--- finished merge after %f seconds\n',mergetime);
                end
                obj.t = iter;
                obj.maxk = iter;
                obj.time = dividetime + conquertime + mergetime;
            else
                pftime = tic;
                set = obj.system.X0;
                obj.sets = HyperRect.empty(1,0);
                % sample initial particles
                X = obj.initialParticles();
                for k = 1:obj.options.iterations
                    z = obj.system.z(:,k+1);
                    u = obj.system.u(:,k);
                    X = obj.runIteration(X,z,u);
                    %set = obj.system.g(set,u) + obj.system.V;
                    %obj.sets(k) = set;
                    if obj.counts(k) == 0
                        break
                    end
                end
                obj.time = toc(pftime);
                obj.cleanup();
            end
            if obj.options.reduceToHull
                obj.reduce();
            end
                
            function updateStatus(data)
                i_ = data(1);
                a_ = data(2);
                p = p + 1;
                if true && a_ > 0
                    fprintf('--- %d/%d: subPF No. %d reached area %f\n',p,n,i_,a_);
                end
            end
        end
        
        function X = initialParticles(obj)
            N = obj.options.samples;
            X = obj.system.X0.sample(N);
        end
        
        function X = runIteration(obj,X,z,u)
            obj.t = obj.t + 1;
            % select algorithm
            switch obj.algorithm
                case PFAlgorithms.Naive
                    X = obj.NaivePF(X,z,u);
                case PFAlgorithms.RedundancyRemoval
                    X = obj.RedundancyRemovalPF(X,z,u);
                case PFAlgorithms.RedundancyRemoval2
                    X = obj.RedundancyRemovalPF2(X,z,u);
                case PFAlgorithms.RedundancyRemoval3
                    X = obj.RedundancyRemovalPF3(X,z,u);
                case PFAlgorithms.WeightOptimisation
                    X = obj.WeightOptimisationPF(X,z,u);
            end
            % remove discarded particles
            obj.particles(:,:,obj.t) = X;
            X = rmmissing(X,2);
            obj.counts(obj.t) = size(X,2);
        end
        
        function X = NaivePF(obj,X,z,u)
        %NAIVEPF Basic Particle Filter Algorithm
        %   
            dim = obj.system.dim;
            N = obj.options.samples;
            N_ = size(X,2);
            % sample propagation
%             X_ = NaN(dim,N_);
%             for i = 1:N_
%                 X_(:,i) = obj.system.g(X(:,i),u);
%             end
            X_ = obj.system.gpredict(X,u);
            % generate new particles
            X = NaN(dim,N);
            for i = 1:N
                i_ = randi([1 N_]); % index selection
                X(:,i) = X_(:,i_) + obj.system.V.sample(); % perturbation
                J = z - obj.system.h(X(:,i)); % innovation calculation
                if ~obj.system.W.contains(J) % discard if not in W
                    X(:,i) = NaN(dim,1);
                end
            end
        end
        
        function X = RedundancyRemovalPF(obj,X,z,u)
        %REDUNDANCYREMOVALPF Particle Filter Algorithm with Redundancy Removal
        %   Use a grid before adding perturbation
            dim = obj.system.dim;
            N = obj.options.samples;
            eps = 0.5*(1/N)^(1/dim); % edge-length of grid cells -> in options
            % sample propagation
            X_ = obj.system.gpredict(X,u);
            % filter particles with grid
            X_ = obj.filterGrid(X_,eps);
            N_ = size(X_,2);
            % generate new particles
            X = NaN(dim,N);
            for i = 1:N
                i_ = randi([1 N_]); % index selection
                X(:,i) = X_(:,i_) + obj.system.V.sample(); % perturbation
                J = z - obj.system.h(X(:,i)); % innovation calculation
                if ~obj.system.W.contains(J) % discard if not in W
                    X(:,i) = NaN(dim,1);
                end
            end
        end
        
        function obj = RedundancyRemovalPF2(obj,X,z,u)
        %REDUNDANCYREMOVALPF Particle Filter Algorithm with Redundancy Removal
        %   Use a grid after adding perturbation
            dim = obj.system.dim;
            N = obj.options.samples;
            N_ = size(X,2);
            eps = (1/N)^(1/dim); % edge-length of grid cells
            % sample propagation
            X_ = NaN(dim,N_);
            for i = 1:N_
                X_(:,i) = obj.system.g(X(:,i),u);
            end
            % generate new particles
            X = NaN(dim,N);
            for i = 1:N
                i_ = randi([1 N_]); % index selection
                X(:,i) = X_(:,i_) + obj.system.V.sample(); % perturbation
                J = z - obj.system.h(X(:,i)); % innovation calculation
                if ~obj.system.W.contains(J) % discard if not in W
                    X(:,i) = NaN(dim,1);
                end
            end
            grid = Grid(X,eps); % setup grid
            % discard particles using grid
            for i = 1:N_
                if grid.freeCell(X(:,i))
                    grid = grid.occupyCell(X(:,i));
                else
                    X(:,i) = NaN(dim,1);
                end
            end
        end
        
        function X = RedundancyRemovalPF3(obj,X,z,u)
        %REDUNDANCYREMOVALPF Particle Filter Algorithm with Redundancy Removal
        %   Use a grid before and after adding perturbation
            dim = obj.system.dim;
            N = obj.options.samples;
            eps = (1/N)^(1/dim); % edge-length of grid cells
            % sample propagation
            N_ = size(X,2);
            X_ = NaN(dim,N_);
            for i = 1:N_
                X_(:,i) = obj.system.g(X(:,i),u);
            end
            grid = Grid(X_,eps); % setup first grid
            % discard particles using first grid
            for i = 1:N_
                if grid.freeCell(X_(:,i))
                    grid = grid.occupyCell(X_(:,i));
                else
                    X_(:,i) = NaN(dim,1);
                end
            end
            % remove discarded particles
            X_ = rmmissing(X_,2);
            N_ = size(X_,2);
            % generate new particles
            X = NaN(dim,N);
            for i = 1:N
                i_ = randi([1 N_]); % index selection
                X(:,i) = X_(:,i_) + obj.system.V.sample(); % perturbation
                J = z - obj.system.h(X(:,i)); % innovation calculation
                if ~obj.system.W.contains(J) % discard if not in W
                    X(:,i) = NaN(dim,1);
                end
            end
            % remove discarded particles
            grid = Grid(X,eps); % setup second grid
            % discard particles using second grid
            for i = 1:N_
                if grid.freeCell(X(:,i))
                    grid = grid.occupyCell(X(:,i));
                else
                    X(:,i) = NaN(dim,1);
                end
            end
        end
        
        function X = filterGrid(~,X,eps)
            % setup grid
            grid = Grid(X,eps);
            % discard particles using grid
            [dim,N] = size(X);
            for i = 1:N
                if grid.freeCell(X(:,i))
                    grid = grid.occupyCell(X(:,i));
                else
                    X(:,i) = NaN(dim,1);
                end
            end
            % remove discarded particles
            X = rmmissing(X,2); 
        end
        
        function X = WeightOptimisationPF(obj,X,z,u)
            dim = obj.system.dim;
            N = obj.options.samples;
            % sample propagation
            N_ = size(X,2);
            X_ = NaN(dim,N_);
            for i = 1:N_
                X_(:,i) = obj.system.g(X(:,i),u);
            end
            if obj.options.debug
                %disp(strcat('k=',num2str(k),': start optimisation'));
            end
            q = obj.optimiseWeights(X_,obj.system.V); % optimise sample weights
            % generate new particles
            X = NaN(dim,N);
            for i = 1:N
                i_ = obj.weightedSelection(1,N_,q); % index selection
                X(:,i) = X_(:,i_) + obj.system.V.sample(); % perturbation
                J = z - obj.system.h(X(:,i)); % innovation calculation
                if ~obj.system.W.contains(J) % discard if not in W
                    X(:,i) = NaN(dim,1);
                end
            end
        end
        
        function q = optimiseWeights(obj,X,V)
        %OPTIMISEWEIGHTS Generates optimised weights given sets X and V
        %   X is the set of particles at current timestamp, V is the support of a rv
            N = size(X,2);
            f = zeros(N+1,1);
            f(N+1) = 1;
            % bounds
            lb = zeros(N+1,1);
            ub = ones(N+1,1);
            % equality constraints
            Aeq = zeros(1,N+1);
            Aeq(1:N) = 1;
            beq = 1;
            % inequality constraints
            A = zeros(N+1);
            A(:,N+1) = -1;
            b = zeros(N+1,1);
            for i = 1:N
                A(i,1:N) = V.contains(X-X(:,i)); % loop over j hidden in vectorization
            end
            %m = floor(N/2)+1;
            %A(m:N,:) = [];
            %b(m:N,:) = [];
            % fprintf('Besetzt mit %f Prozent\n',sum(A(:) > 0)/numel(A)*100);
            %imagesc(A);
            % solve linear program
            A = sparse(A);
            if obj.options.gurobi
                x = gurobi_linprog(f,A,b,Aeq,beq,lb,ub,obj.options.LPOptions);
            else
                x = linprog(f,A,b,Aeq,beq,lb,ub,obj.options.LPOptions);
            end
            if numel(x) == N+1
                %y = x(N+1);
                q = x(1:N)/sum(x(1:N));
            else
                q = ones(N,1)/N; % uniform weights when no solution found
            end
        end
        
        function i = weightedSelection(~,a,b,q)
        %WEIGHTEDSELECTION Weighted selection of an integer from the interval [a,b]
        %   Uses the weights in the array q
            r = rand(1,1); % generate random number r
            S = cumsum(q);
            ind = a:b; % available indices
            i = min(ind(S>r)); % choose the index according to r
            if isempty(i)
                i = b;
            end
        end
        
        function plot(obj, showParticles, showConvhull, showSet, autoScale, fig)
            %PLOT Plot the particles with convex hull
            %   if showParticles is false: only convex hull
            if nargin < 6
                figure
                if nargin < 5
                    autoScale = true;
                    if nargin < 4
                        showSet = false;
                        if nargin < 3
                            showConvhull = true;
                            if nargin < 2
                                showParticles = true;
                            end
                        end
                    end
                end
            else
                figure(fig);
            end
            set(gcf, 'Position', get(0, 'Screensize'));
            % plot set of feasible states
            limits = obj.particleLimits();
            rows = ceil(sqrt(max(obj.t,obj.options.iterations)));
            for k = 1:obj.t
                subplot(rows,rows,k);
                hold on;
                if showSet
                    obj.plotIntv(k);
                end
                if showParticles
                    obj.plotParticles(k);
                end
                obj.plotX(k);
%                 if showSet && numel(obj.sets) >= k
%                     obj.sets(k).plot();
%                 end
                if showConvhull
                    area = obj.plotConvhull(k);
                    title(strcat('k=',num2str(k),', area=',num2str(area)));
                else
                    title(strcat('k=',num2str(k)));
                end
                if autoScale
                    xlim(limits(1,:))
                    ylim(limits(2,:))
                end
                hold off;
            end
        end
        
        function plot2(obj, showParticles, showConvhull, showSet, autoScale, fig)
            %PLOT2 Plot the particles with convex hull
            %   if showParticles is false: only convex hull
            if nargin < 6
                figure
                if nargin < 5
                    autoScale = true;
                    if nargin < 4
                        showSet = false;
                        if nargin < 3
                            showConvhull = true;
                            if nargin < 2
                                showParticles = true;
                            end
                        end
                    end
                end
            else
                figure(fig);
            end
            set(gcf, 'Position', get(0, 'Screensize'));
            % plot set of feasible states
            limits = obj.particleLimits();
            rows = ceil(sqrt(max(obj.t,obj.options.iterations)));
            for k = 1:obj.t
                subplot(rows,rows,k);
                hold on;
                if showSet
                    obj.plotIntv2(k);
                end
                if showParticles
                    obj.plotParticles2(k);
                end
                obj.plotX2(k);
%                 if showSet && numel(obj.sets) >= k
%                     obj.sets(k).plot();
%                 end
                if showConvhull
                    area = obj.plotConvhull2(k);
                    title(strcat('k=',num2str(k),', area=',num2str(area)));
                else
                    title(strcat('k=',num2str(k)));
                end
                if autoScale
                    xlim(limits(3,:))
                    ylim(limits(4,:))
                end
                hold off;
            end
        end
        
        function plotParticles(obj,k)
            P = obj.particlesAt(k);
            scatter(P(1,:),P(2,:),3,'MarkerEdgeColor','b','MarkerEdgeAlpha',.2);
        end
        
        function plotParticles2(obj,k)
            P = obj.particlesAt(k);
            scatter(P(3,:),P(4,:),3,'MarkerEdgeColor','b','MarkerEdgeAlpha',.2);
        end
        
        function crossplotParticles(obj,k,dim,showSets,autoScale)
            if nargin < 5
                autoScale = true;
                if nargin < 4
                    showSets = false;
                end
            end
            figure;
            sgtitle(['cross plot particles at timestep k = ',num2str(k)]);
            P = obj.particlesAt(k);
            m = size(P,1);
            if nargin == 3
                m = min(m,dim);
            end
            limits = obj.particleLimits(k,0.01);
            for i = 1:m
                for j = 1:i
                    subplot(m,m,j+(i-1)*m);
                    if i == j
                        histogram(P(i,:),'FaceColor','b');
                        title(['histogram x',num2str(i)]);
                        xlabel(['x',num2str(i)]);
                    else
                        hold on
                        if showSets
                            if isa(obj.intvE,'cell')
                                E = obj.intvE{k+1};
                                if size(E,1) > 0
                                    plotIntervals(E(:,[j i],:),'g',0);
                                end
                            end
                            if isa(obj.intvS,'cell')
                                S = obj.intvS{k+1};
                                if size(S,1) > 0
                                    plotIntervals(S(:,[j i],:),'g',0);
                                end
                            end
                            if isa(obj.intvSz,'cell')
                                Sz = obj.intvSz{k+1};
                                if size(Sz,1) > 0
                                    plotIntervals(Sz,'c',0);
                                end
                            end
                        end
                        scatter(P(j,:),P(i,:),3,'MarkerEdgeColor','b','MarkerEdgeAlpha',.15);
                        x = obj.system.x(:,k+1); %k+1 for MPC!
                        scatter(x(j),x(i),6,'MarkerEdgeColor','r','MarkerFaceColor','r');
                        hold off
                        if autoScale
                            xlim(limits(j,:));
                            ylim(limits(i,:));
                        end
                        title(['x',num2str(j),' & x',num2str(i)]);
                        xlabel(['x',num2str(j)]);
                        ylabel(['x',num2str(i)]);
                    end
                end
            end
        end
        
        function area = plotConvhull(obj,k)
            if obj.counts(k) == 0
                area = 0;
            else
                [c,area] = obj.convhull(k);
                plot(c(1,:),c(2,:),'k');
            end
        end
        
        function area = plotConvhull2(obj,k)
            if obj.counts(k) == 0
                area = 0;
            else
                [c,area] = obj.convhull2(k);
                plot(c(1,:),c(2,:),'k');
            end
        end
        
        function area = plotBoundary(obj,k)
            if obj.counts(k) == 0
                area = 0;
            else
                [c,area] = obj.boundary(k);
                plot(c(1,:),c(2,:),'k');
            end
        end
        
        function area = plotBoundary2(obj,k)
            if obj.counts(k) == 0
                area = 0;
            else
                [c,area] = obj.boundary(k);
                plot(c(1,:),c(2,:),'k');
            end
        end
        
        function x = plotX(obj,k)
            x = obj.system.x(:,k+1); %k+1 for MPC!
            scatter(x(1),x(2),6,'MarkerEdgeColor','r','MarkerFaceColor','r');
        end
        
        function x = plotX2(obj,k)
            x = obj.system.x(:,k+1); %k+1 for MPC!
            scatter(x(3),x(4),6,'MarkerEdgeColor','r','MarkerFaceColor','r');
        end
        
        function plotXRange(obj,k)
            x = obj.system.x(:,k+1);
            w = obj.system.W.bounds;
            rectangle('Position',[x(1)+w(1,1),x(2)+w(2,1),w(1,2)-w(1,1),w(2,2)-w(2,1)],'EdgeColor','red');
        end
        
        function x = plotZ(obj,k)
            x = obj.system.z(:,k+1);
            scatter(x(1),x(2),6,'MarkerEdgeColor','g','MarkerFaceColor','g');
        end
        
        function plotZRange(obj,k)
            z = obj.system.z(:,k+1);
            w = obj.system.W.bounds;
            rectangle('Position',[z(1)+w(1,1),z(2)+w(2,1),w(1,2)-w(1,1),w(2,2)-w(2,1)],'EdgeColor','green');
        end
        
        function plotIntv(obj,k)
            if isa(obj.intvE,'cell')
                E = obj.intvE{k+1};
                if size(E,1) > 0
                    plotIntervals(E(:,1:2,:),'g',0);
                end
            end
            if isa(obj.intvS,'cell')
                S = obj.intvS{k+1};
                if size(S,1) > 0
                    plotIntervals(S(:,1:2,:),'g',0);
                end
            end
            if isa(obj.intvSz,'cell')
                Sz = obj.intvSz{k+1};
                if size(Sz,1) > 0
                    plotIntervals(Sz,'c',0);
                end
            end
%             obj.drawIntv(S,[1 0 0 0.1]);
%             obj.drawIntv(E,[1 0 0 1]);
        end
        
        function plotIntv2(obj,k)
            S = obj.intvS{k+1};
            E = obj.intvE{k+1};
            if size(E,1) > 0
                plotIntervals(E(:,3:4,:),'y',0);
            end
            if size(S,1) > 0
                plotIntervals(S(:,3:4,:),'g',0);
            end
%             obj.drawIntv(S,[1 0 0 0.1]);
%             obj.drawIntv(E,[1 0 0 1]);
        end
        
        function drawIntv(obj,I,color)
            x = I(:,1);
            y = I(:,2);
            arrayfun(@(xl,xu,yl,yu) rectangle('Position',[xl yl xu-xl yu-yl],'LineWidth',0.5,'LineStyle','-','EdgeColor',color,'FaceColor',color), x.lower, x.upper, y.lower, y.upper);
        end
        
        function plotRanges(obj,k)
            hold on
            obj.plotX(k);
            obj.plotZ(k);
            obj.plotXRange(k);
            obj.plotZRange(k);
            obj.plotConvhull(k);
        end
        
        function r = particleLimits(obj,k,eps)
            if nargin < 3
                eps = 0;
                if nargin < 2
                    k = 1:max(obj.maxk,obj.t);
                end
            end
            P = obj.particles(:,:,k);
            r = [min(P,[],[2 3])-eps max(P,[],[2 3])+eps];
        end
        
        function v = volume(obj,k)
        %CONVHULL Get convex hull around particles
		%   for iteration/timestep k
            try
                P = obj.particlesAt(k);
                [c_,v_] = convhulln(P');
                %c = P(:,c_);
                v = v_;
            catch
                %c = [];
                v = 0;
            end
        end
        
        function v = boxVolume(obj,k)
            try
                P = obj.particlesAt(k);
                v = prod(max(P,[],2)-min(P,[],2));
            catch
                v = 0;
            end
        end
        
        function [c,a] = convhull(obj,k)
        %CONVHULL Get convex hull around particles
		%   for iteration/timestep k
            try
                P = obj.particlesAt(k);
                P = P(1:2,:); % only 2D!
                [c_,a_] = convhull(P');
                c = P(:,c_);
                a = a_;
            catch
                c = [];
                a = 0;
            end
        end
        
        function [c,a] = convhull2(obj,k)
        %CONVHULL Get convex hull around particles
		%   for iteration/timestep k
            try
                P = obj.particlesAt(k);
                P = P(3:4,:); % only 2D!
                [c_,a_] = convhull(P');
                c = P(:,c_);
                a = a_;
            catch
                c = [];
                a = 0;
            end
        end
        
        function [c,a] = boundary(obj,k)
        %BOUNDARY Get boundary around particles
		%   for iteration/timestep k
            try
                P = obj.particlesAt(k);
                P = P(1:2,:); % only 2D!
                [c_,a_] = boundary(P');
                c = P(:,c_);
                a = a_;
            catch
                c = [];
                a = 0;
            end
        end
        
        function a = area(obj,k)
        %AREA Get area of convex hull around particles
		%   for iteration/timestep k
            [~,a_] = obj.convhull(k);
            a = a_;
        end
        
        function reduce(obj)
        %REDUCE Reduce particles to their convex hull
		%   can be called to free up memory space
            [d,s,iter] = size(obj.particles);
            for k = 1:iter
                P = obj.convhull(k);
                r = size(P,2);
                obj.counts(k) = r;
                obj.particles(:,1:r,k) = P;
                obj.particles(:,r+1:s,k) = NaN(d,s-r);
            end
            obj.cleanup();
        end
        
        function P = particlesAt(obj,k,n)
        %PARTICLESAT Get n (default: all) particles
		%   for iteration/timestep k
            P = rmmissing(obj.particles(:,:,k),2);
            if nargin == 3
                [dim,N] = size(P);
                if n < N
                    eps = (obj.boxVolume(k)/n)^(1/dim);
                    P = obj.filterGrid(P,eps);
                    N = size(P,2);
                end
                P = P(:,randperm(N,min(n,N)));
            end
        end
        
        function range = particleRangeAt(obj,k)
            if nargin < 2
                k = 0:obj.t;
            end
            if numel(k) > 1
                range = [];
                for i = 1:numel(k)
                    range = cat(3,range,obj.particleRangeAt(k(i)));
                end
            else
                if k == 0
                    range = obj.system.X0.bounds;
                else
                    P = obj.particlesAt(k);
                    range = [min(P,[],2) max(P,[],2)];
                end
            end
        end
        
        function merge(obj,other)
        %MERGE Merge particles of other into Particle Filter
        %   other is also a object of class ParticleFilter
            if numel(other) > 1
                for i = 1:numel(other)
                    obj.merge(other(i));
                end
            else
                if max(other.counts) > 0
                    [d1,s1,k1] = size(obj.particles);
                    [d2,s2,k2] = size(other.particles);
                    if k1 ~= k2 % pad if not same size
                        p1 = max([k1,k2]) - k1;
                        p2 = max([k1,k2]) - k2;
                        obj.particles = cat(3,obj.particles,NaN(d1,s1,p1));
                        other.particles = cat(3,other.particles,NaN(d2,s2,p2));
                    end
                    obj.particles = cat(2,obj.particles,other.particles);
                    obj.counts = obj.counts + other.counts;
                end
            end
        end
        
        function clear(obj)
        %CLEAR Clear particles etc.
        %   should be called to free memory space
            obj.particles = [];
            obj.counts = [];
        end
        
        function cleanup(obj)
        %CLEANUP Free up unused memory space after `run`
        %   should be called to free memory space
            obj.maxk = size(obj.particles,3);
            maxCount = max(obj.counts) + 1;
            if isnan(maxCount)
                maxCount = 1;
            end
            obj.particles(:,maxCount:size(obj.particles,2),:) = [];
        end
        
        function p = hasParallel(~)
            try
                p = license('test','Distrib_Computing_Toolbox');
            catch
                p = false;
            end
        end
        
    end
end

