classdef MPC < handle
    %MPC Robust Model Predictive Control
    %   based on a Set-Membership Particle Filter
    
    properties
        system % system to be controlled
        mode % MPC mode (0: measurement, 1: PF estimation, 2: exact state)
        PF % particle filter
        PFdraws % number of particles drawn from PF
        H % prediction horizon
        Q % state penalty matrix
        R % input penalty matrix
        P % final state penalty matrix
        u_lb % control input lower bound
        u_ub % control input upper bound
        r % reference function
        r_vals % reference values
        x % real system state
        z % system measurements
        X % estimated system states
        u % calculated control inputs
        u_last % last control input
        estX_last % last state estimation
        v % process noise
        w % measurement noise
        noiseScale % noise scale factor
        simulate % use random noise (logical)
        convhull % use only convhull of particles (logical)
    end
    
    methods
        function obj = MPC(system,mode,samples)
            %MPC Construct an instance of this class
            %   Initialize system and Particle Filter
            if nargin < 3
                samples = 50000;
                if nargin < 2
                    mode = 0;
                end
            end
            obj.system = system;
            obj.Q = diag([1 1 0.01 0.01]);
            obj.R = 1e-5 * eye(obj.system.dimu);
            obj.P = zeros(obj.system.dim);
            obj.u_lb = [-10;-10];
            obj.u_ub = [10;10];
            logRef = @(t) [pi;pi;0;0] .* 1./(1+exp(-0.2*t+5)) + [0;0;pi;pi] .* 0.2*exp(-0.2*t+5)./(1+exp(-0.2*t+5)).^2;
            obj.r = @(t) [pi;pi;0;0] + 0 * t; % reference function
%             obj.r = @(t) logRef(t);
%             obj.r = @(t) [pi;pi;0;0] + pi/8 * [sin(2*pi*t/40);sin(2*pi*t/40);40/(2*pi)*cos(2*pi*t/40);40/(2*pi)*cos(2*pi*t/40)];
%             obj.r = @(t) wrapToPi(system.x0 + obj.system.T*[pi/8;pi/8;0;0].*t); %[pi;pi;0;0] + [-pi/4;pi/4;0;0] * sin(2*pi*t/40);
            obj.x = system.x0;
            obj.u = [];
            obj.v = [];
            obj.w = [];
            obj.H = 10;
            obj.noiseScale = 1.0;
            obj.mode = mode;
            obj.simulate = true;
            obj.convhull = false;
            obj.PFdraws = 1000;
            % initialize Particle Filter
            PFoptions = PFOptions();
            PFoptions.samples = samples;
            obj.PF = ParticleFilter(obj.system,PFAlgorithms.Naive,PFoptions);
        end
        
        function c = cost(obj,x0,u,xref)
            %COST Evaluate MPC cost functional for all x0
            %   x0: array of start system states, u: control inputs
            if isa(x0,'interval')
                %c = obj.costIntv(varargin{:});
            else
                dim = obj.system.dim;
                if min(size(u) ~= [dim obj.H]) 
                    error('wrong dimensions');
                end
                if nargin < 4
                    xref = zeros(obj.system.dim,obj.H+1);
                end
                % predict system states for all x0
                xt = obj.system.predict(x0,u);
                % evaluate MPC cost functional
                uref = obj.system.uref(xref);
                c = obj.Jv(xt,u,xref,uref);
                % choose worst possible system state
                c = max(c);
            end
        end
        
        function c = costIntv(obj,x0,u,xref)
            %COST Evaluate MPC cost functional for all x0
            %   x0: array of start system states, u: control inputs
            dim = obj.system.dim;
            if min(size(u) ~= [dim obj.H]) 
                error('wrong dimensions');
            end
            if nargin < 4
                xref = zeros(obj.system.dim,obj.H+1);
            end
            % predict system states for all x0
            xt = obj.system.predictIntv(x0,u);
            % evaluate MPC cost functional
            uref = obj.system.uref(xref);
            c = obj.JIntv(xt,u,xref,uref);
            % choose worst possible system state
            c = max(c.upper(:));
        end
        
        function c = J(obj,xt,u)
            %J MPC cost function
            %   xt: array of system states, u: control inputs
            N = size(xt,3);
            c = zeros(1,N);
            for i = 1:N
                for k = 1:obj.H
                    c(i) = c(i) + xt(:,k,i)' * obj.Q * xt(:,k,i) + u(:,k)' * obj.R * u(:,k);
                end
                c(i) = c(i) + xt(:,obj.H+1,i)' * obj.P * xt(:,obj.H+1,i);
            end
        end
        
        function c = Jv(obj,xt,u,xref,uref)
            %Jv Vectorized MPC cost function
            %   xt: array of system states, u: control inputs
            if nargin < 4
                xref = zeros(obj.system.dim,obj.H+1);
            end
            N = size(xt,3);
            c = zeros(N,1);
            for k = 1:obj.H
                %blkdiag
                xk = squeeze(permute(xt(:,k,:),[3 1 2])) - xref(:,k)';
                uk = u(:,k) - uref(:,k);
                % uref bestimmen, Differenz u - uref bilden
                c = c + sum(xk * obj.Q .* xk,2) + uk' * obj.R * uk;
            end
            xH = squeeze(permute(xt(:,obj.H+1,:),[3 1 2])) - xref(:,obj.H+1)';
            c = c + sum(xH * obj.P .* xH,2);
        end
        
        function c = JIntv(obj,xt,u,xref,uref)
            %JIntv MPC cost function for intervals
            %   xt: intervals which contain system states, u: control inputs
            if nargin < 4
                xref = zeros(obj.system.dim,obj.H+1);
            end
            N = size(xt,3);
            c = zeros(N,1);
            for k = 1:obj.H
                %blkdiag
                xk = (xt{k}' - xref(:,k))';
                uk = u(:,k) - uref(:,k);
                % uref bestimmen, Differenz u - uref bilden
                c = c + quadr(obj.Q,xk) + uk' * obj.R * uk;
            end
            xH = (xt{obj.H+1}' - xref(:,obj.H+1))';
            c = c + quadr(obj.P,xH);
        end
        
        function c = realJ(obj,normalize)
            if nargin < 2
                normalize = false;
            end
            c = 0;
            for k = 1:size(obj.x,2)-1
                xk = obj.r_vals(:,k) - obj.x(:,k);
                c = c + xk' * obj.Q * xk + obj.u(:,k)' * obj.R * obj.u(:,k);
            end
            xH = obj.r_vals(:,size(obj.x,2)) - obj.x(:,size(obj.x,2));
            c = c + xH' * obj.P * xH;
            if normalize
                c = c * obj.H / size(obj.x,2);
            end
        end
        
        function u = run(obj,t)
			%RUN Run MPC for timestep t
			%   
            if t == 0
                obj.x(:,1) = obj.system.x0;
            end
            
            % measurement
            if obj.simulate
                obj.w(:,t+1) = obj.noiseScale .* obj.system.Wsim.sample();
            elseif size(obj.w,2) < t+1
                obj.w(:,t+1) = zeros(1,obj.system.dimz);
            end
            obj.z(:,t+1) = obj.system.h(obj.x(:,t+1)) + obj.w(:,t+1);
            
            % state estimation by Particle Filter
            estX = obj.estimateState(t);
            
            % state reference
            ref = obj.r(t:t+obj.H);
            
            % find u by optimization
            if t == 0
                u0 = zeros(obj.system.dimu,obj.H);
            else
                u0 = [obj.u_last(:,2:end), obj.u_last(:,end)];
            end
            options = optimoptions('fmincon','Display','off','Algorithm','sqp');
            problem.options = options;
            problem.solver = 'fmincon';
            problem.objective = @(u)obj.cost(estX, u, ref);
            problem.x0 = u0;
            problem.lb = obj.u_lb .* ones(1,obj.H);
            problem.ub = obj.u_ub .* ones(1,obj.H);
            u = fmincon(problem);
            %u = u0;
            
            % apply u to system
            if obj.simulate
                obj.v(:,t+1) = obj.noiseScale .* obj.system.Vsim.sample();
            elseif size(obj.v,2) < t+1
                obj.v(:,t+1) = zeros(1,obj.system.dim);
            end
            obj.u(:,t+1) = u(:,1);
            obj.u_last = u;
            obj.x(:,t+2) = obj.system.g(obj.x(:,t+1),obj.u(:,t+1)) + obj.v(:,t+1);
            obj.r_vals(:,t+1:t+obj.H+1) = ref;
            obj.system.x = obj.x;
            obj.system.z = obj.z;
            obj.PF.system = obj.system;
            J = obj.realJ(true);
            disp(['t = ',num2str(t),', J = ',num2str(J)]);
        end
        
        function u = runUntil(obj,Tend,Tstart)
            if nargin < 3
                Tstart = 0;
            end
            for t = Tstart:Tend
                obj.run(t);
            end
            u = obj.u;
        end
        
        function estX = estimateState(obj,t)
            switch obj.mode
                case 0
                    zt = obj.z(:,t+1);
%                     if t == 0
%                         estX = zeros(obj.system.dim,1);
%                     else
%                         estX = obj.estX_last;
%                     end
%                     fun = @(x)obj.system.h(x) - zt;
%                     %xt = fsolve(fun,estX);
%                     options = optimoptions('fmincon','Display','off');
%                     problem.options = options;
%                     problem.solver = 'fmincon';
%                     problem.objective = @(x)sum(abs(fun(x)));
%                     problem.x0 = estX;
%                     problem.lb = estX - [-2*pi;-2*pi;-100;-100];
%                     problem.ub = [100;100;100;100];
%                     xt = fmincon(problem);
%                     estX(1:numel(xt)) = xt;
                    estX = zt;
                case 1
                    if t == 0
                        estX = obj.system.x0;
                        obj.X = obj.PF.initialParticles();
                    else
                        obj.X = obj.PF.runIteration(obj.X, obj.z(:,t+1), obj.u(:,t));
                        if obj.convhull
                            estX = obj.PF.convhull(t);
                        else
                            estX = obj.PF.particlesAt(t,obj.PFdraws);
                        end
                        if numel(estX) == 0
                            error('no particles at timestep %t\n', t);
                        end
                    end
                case 2
                    xt = obj.x(:,t+1);
                    estX = zeros(obj.system.dim,1);
                    estX(1:numel(xt)) = xt;
            end
            obj.estX_last = estX;
        end
        
        function plotStates(obj, fig)
            if nargin < 2
                fig = 1;
            end
            figure(fig);
            title('Systemzustände');
            hold on
            plot(obj.x(1,1:end-1));
            plot(obj.x(2,1:end-1));
            legend('x1','x2');
            xlim([1 200]);
            hold off
        end
        
        function plot(obj, showRef, showMeas, showUncty, fig)
            %PLOT Plot the system states and control inputs
            %   fig: figure number
            if nargin < 5
                figure();
                if nargin < 4
                    showUncty = false;
                    if nargin < 3
                        showMeas = false;
                        if nargin < 2
                            showRef = true;
                        end
                    end
                end
            else
                figure(fig);
            end
            %set(gcf, 'Position', get(0, 'Screensize'));
            % plot system states
            subplot(2,1,1);
            hold on
            title('Systemzustände');
            t = 0:(size(obj.x,2)-2);
%             yyaxis left
            W = obj.system.W;
            if showUncty
                switch obj.mode
                    case 0
                        obj.plotShade(t,obj.x(1,1:end-1)+W.bounds(1,1),obj.x(1,1:end-1)+W.bounds(1,2),'b',0.05);
                        obj.plotShade(t,obj.x(2,1:end-1)+W.bounds(2,1),obj.x(2,1:end-1)+W.bounds(2,2),'r',0.05);
                    case 1
                        obj.plotShade(t,obj.x(1,1:end-1)+W.bounds(1,1),obj.x(1,1:end-1)+W.bounds(1,2),'b',0.05);
                        obj.plotShade(t,obj.x(2,1:end-1)+W.bounds(2,1),obj.x(2,1:end-1)+W.bounds(2,2),'r',0.05);
                        PRange = obj.PF.particleRangeAt();
                        obj.plotShade(t,PRange(1,1,:),PRange(1,2,:),'b',0.2);
                        obj.plotShade(t,PRange(2,1,:),PRange(2,2,:),'r',0.2);
                end
            end
            if showMeas
                scatter(t,obj.z(1,1:numel(t)),10,'b','+','HandleVisibility','off');
                scatter(t,obj.z(2,1:numel(t)),10,'r','+','HandleVisibility','off');
            end
            plot(t,obj.x(1,1:end-1),'b-');
            plot(t,obj.x(2,1:end-1),'r-');
            legend({'$x_1$','$x_2$'},'Interpreter','latex','NumColumns',2);
%             yyaxis right
%             plot(t,obj.x(3,1:end-1),'b-.');
%             plot(t,obj.x(4,1:end-1),'r-.');
%             legend('x1','x2','x3','x4');
            % reference values
            if showRef
                plot(t,obj.r_vals(1,1:size(obj.x,2)-1),'b--');
                plot(t,obj.r_vals(2,1:size(obj.x,2)-1),'r--');
            else
                yline(pi,'-.k','$x_\mathrm{ref}$','Interpreter','latex','HandleVisibility','off');
            end
            %ylim([2.5 3.5]);
            %xlim([0 6000]);
            hold off
            % plot control inputs
            subplot(2,1,2);
            hold on
%             title('total energy');
            title('Stellgrößen');
%             plot(obj.system.energy(obj.x))
            stairs(t,obj.u(1,:),'b-');
            stairs(t,obj.u(2,:),'r-');
            legend({'$u_1$','$u_2$'},'Interpreter','latex','NumColumns',2);
            ylim([min(obj.u_lb) max(obj.u_ub)]);
            %xlim([0 6000]);
            hold off
        end
        
        function plotShade(obj,t,x1,x2,c,a)
            S = [t, fliplr(t)];
            inBetween = [reshape(x1,1,[]), fliplr(reshape(x2,1,[]))];
            fill(S, inBetween, c,'FaceColor', c, 'FaceAlpha', a, 'EdgeColor','none','HandleVisibility','off');
        end
    end
end

