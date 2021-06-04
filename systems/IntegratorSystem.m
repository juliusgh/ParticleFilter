classdef IntegratorSystem < System
    %SYSTEM Summary of this class goes here
    %   Detailed explanation goes here
    methods
        function obj = IntegratorSystem()
            %SYSTEM Construct an instance of this class
            %   Simple mechanical test system
            obj.dim = 2; % system dimensions
            obj.dimz = 2; % measurement dimensions
            obj.dimu = 2; % control input dimensions
            obj.x0 = [1;0]; % initial condition
            obj.X0 = obj.x0 + HyperRect([-0.5 0.5],[-0.5 0.5]); % initial noise support
            obj.V = HyperRect([-0.03 0.03],[-0.03 0.03]); % system noise support
            obj.W = HyperRect([-0.3 0.3],[-0.5 0.5]); % measurement noise support
            obj.uk = @(k) [0;0]; % fixed input function
        end
        
        function x_ = g(~,x,u)
        %G process function with input
        %   x: current system state, u: control input
            x_ = [1,1;0,1]*x + [0,0;0,1]*u;
        end
        
        function z = h(~,x)
        %G measurement function
        %   x: current system state
            z = x;
        end
    end
end

