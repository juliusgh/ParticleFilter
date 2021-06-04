classdef System
    %SYSTEM Base class for dynamical systems
    %   Contains test system for Particle Filters
    
    properties
        dim % system dimensions
        dimz % measurement dimensions
        dimu % control input dimensions
        x0 % initial system state
        X0 % initial noise support
        V % system noise support
        W % measurement noise support
        Vsim % system noise support
        Wsim % measurement noise support
        uk % fixed input function
        x % system states
        z % measurements
        u % inputs
    end
    methods
        function obj = System()
            %SYSTEM Construct an instance of this class
            %   Default test system for Particle Filters
            obj.dim = 2; % system dimensions
            obj.dimz = 1; % measurement dimensions
            obj.dimu = 2; % control input dimensions
            obj.x0 = [3;1]; % initial condition
            obj.X0 = HyperRect([0 5],[0 5]); % initial noise support
            obj.V = HyperRect([-0.03 0.03],[-0.03 0.03]); % system noise support
            obj.W = HyperRect([-0.3 0.3]); % measurement noise support
            obj.Vsim = HyperRect([-0.03 0.03],[-0.03 0.03]); % system noise support
            obj.Wsim = HyperRect([-0.3 0.3]); % measurement noise support
            obj.uk = @(k) [0;0] * sin(2*pi*k/8); % fixed input function
        end
        
        function x_ = g(~,x,u)
        %G process function with input
        %   x: current system state, u: control input
            if nargin < 3
                u = [0;0];
            end
            if isa(x,'interval')
                x1 = x(:,1);
                x_ = x + u;
                k = 0.16 .* 0.1 .* (x1.^2);
                x_(:,1) = x_(:,1) + k .* -2;
                x_(:,2) = x_(:,2) + k .* 1;
            else
                x1 = x(1,:,:);
                x_ = x + 0.16 * 0.1 * (x1.^2) .* [-2;1] + u;
            end
        end
        
        function x_ = gpredict(obj,x,u)
        %G process function with input for prediction
        %   x: current system state, u: control input
            x_ = obj.g(x,u);
        end
        
        function z = h(~,x)
        %G measurement function
        %   x: current system state
            if isa(x,'interval')
                z = x*[1,1];
            else
                z = [1,1]*x;
            end
        end
        
        function u = uref(obj,xref)
            %UREF get desired control input for reference state
            %   x: reference system state
            xr = obj.f(0,xref,[0;0]);
            u = -xr;
        end
        
        function obj = simulate(obj,iterations,noise,gfun)
            %SIMULATE Simulate system states and measurements
            %   with perturbation sampled from V and W
            if nargin < 4
                gfun = @obj.g;
                if nargin < 3
                    noise = true;
                end
            end
            for k = 1:iterations+1
                v = 0;
                w = 0;
                if noise
                    v = obj.Vsim.sample();
                    w = obj.Wsim.sample();
                end
                obj.u(:,k) = obj.uk(k-1); % calculate input
                if k == 1
                    obj.x(:,k) = obj.x0;
                else
                    obj.x(:,k) = gfun(obj.x(:,k-1),obj.u(:,k)) + v; % calculate system state
                end
                obj.z(:,k) = obj.h(obj.x(:,k)) + w; % calculate measurement
            end
        end
        
        function obj = simulatePrediction(obj,iterations,noise)
            %SIMULATE Simulate system states and measurements
            %   with perturbation sampled from V and W
            if nargin < 3
                noise = true;
            end
            obj = obj.simulate(iterations,noise,@obj.gpredict);
        end
        
        function xt = predict(obj,x0,u)
		%PREDICT Predict the prospective system states given 
		%   x0: initial states, u: control inputs
            if size(u,1) ~= obj.dimu || size(x0,1) ~= obj.dim
                error('wrong dimensions');
            end
            H = size(u,2); % prediction horizon
            N = size(x0,2); % count of initial states
            xt = NaN(obj.dim,N,H+1);
            xt(:,:,1) = x0;
            for k = 1:H
                xt(:,:,k+1) = obj.gpredict(xt(:,:,k),u(:,k)); % propagate states
            end
            xt = permute(xt, [1 3 2]); % permute for more convenient indexing
        end
        
        function xt = predictIntv(obj,X0,u)
		%PREDICT Predict the prospective system states given 
		%   x0: initial states, u: control inputs
            if size(u,1) ~= obj.dimu || size(X0,2) ~= obj.dim
                error('wrong dimensions');
            end
            H = size(u,2); % prediction horizon
            xt = cell(1,H+1);
            xt{1} = X0;
            for k = 1:H
                xt{k+1} = obj.gpredict(xt{k},u(:,k)); % propagate states
            end
            %xt = permute(xt, [1 3 2]); % permute for more convenient indexing
        end
        
        function parts = divide(obj,divisions)
            %DIVIDE Divide system along every axis in `divisions` parts
            %   for Particle Filters using the divide and conquer principle
            subX0 = obj.X0.divide(divisions);
            n = numel(subX0);
%             parts = System.empty(0,n);
            for i = 1:n
                parts(i) = obj;
                parts(i).X0 = subX0(i);
            end
        end
        
        function plot(obj, fig)
            %PLOT Plot the system states
            %   fig: figure number
            if nargin == 2
                figure(fig);
            else
                figure;
            end
            %set(gcf, 'Position', get(0, 'Screensize'));
            % plot system states
            hold on
            title('system states');
            t = 0:(size(obj.x,2)-1);
            plot(t,obj.x(1,1:end),'b-');
            plot(t,obj.x(2,1:end),'r-');
            hold off
        end
    end
end

