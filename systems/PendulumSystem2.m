classdef PendulumSystem2 < PendulumSystem
    %PENDULUMSYSTEM Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function obj = PendulumSystem2()
            %SYSTEM Construct an instance of this class
            %   Simple mechanical test system
            obj.dim = 4; % system dimensions
            obj.dimz = 4; % measurement dimensions
            obj.dimu = 2; % control input dimensions
            obj.x0 = [pi/2;0;0;0]; % initial condition
            obj.V = HyperRect([-0.001 0.001],[-0.001 0.001],[-0.01 0.01],[-0.01 0.01]); % system noise support
            obj.W = HyperRect([-0.3 0.3],[-0.3 0.3],[-0.3 0.3],[-0.3 0.3]); % measurement noise support
            obj.Vsim = HyperRect([-0.001 0.001],[-0.001 0.001],[-0.01 0.01],[-0.01 0.01]); % system noise support
            obj.Wsim = HyperRect([-0.3 0.3],[-0.3 0.3],[-0.3 0.3],[-0.3 0.3]); % measurement noise support
            obj.X0 = obj.x0 + obj.W; % initial noise support
            obj.uk = @(k) [0;0]; % fixed input function
            obj.T = 0.01;
            obj.G = 9.81;
            obj.L = 1;
            obj.m = 1;
            obj.J = obj.m * obj.L^2;
        end
        
        function z = h(obj,x)
            %H measurement function
            %   x: current system state
            [x1,x2,x3,x4] = unpack(x);
            z = [obj.L .* sin(x1) + obj.L .* sin(x2);...
                -obj.L .* cos(x1) - obj.L .* cos(x2);...
                x3 .* cos(x1) + x4 .* cos(x2);...
                x3 .* sin(x1) + x4 .* sin(x2)];
%             z = [obj.L * sin(x1),...
%                 -obj.L * cos(x1),...
%                 obj.L * sin(x1) + obj.L * sin(x2),...
%                 -obj.L * cos(x1) - obj.L * cos(x2)];
%             if ~isa(x,'interval')
%                 z = z';
%             end
        end
    end
end

