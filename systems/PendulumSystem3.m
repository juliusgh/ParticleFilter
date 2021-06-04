classdef PendulumSystem3 < PendulumSystem
    %PENDULUMSYSTEM Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function obj = PendulumSystem3()
            %SYSTEM Construct an instance of this class
            %   Simple mechanical test system
            obj.dim = 4; % system dimensions
            obj.dimz = 4; % measurement dimensions
            obj.dimu = 2; % control input dimensions
            obj.x0 = [pi/2;pi;0;0]; % initial condition
            obj.V = HyperRect([-0.001 0.001],[-0.001 0.001],[-0.011 0.011],[-0.011 0.011]); % system noise support
%             obj.W = HyperRect([-0.1 0.1],[-0.1 0.1],[-0.1 0.1],[-0.1 0.1]); % measurement noise support
            obj.Vsim = HyperRect([-0.0 0.0],[-0.0 0.0],[-0.01 0.01],[-0.01 0.01]); % system noise support
            obj.W = HyperRect([-0.05 0.05],[-0.05 0.05],[-0.05 0.05],[-0.05 0.05]); % measurement noise support
            obj.Wsim = obj.W;
            %obj.Vsim = HyperRect([-0.0 0.0],[-0.0 0.0],[-0.0 0.0],[-0.0 0.0]); % system noise support
            obj.X0 = obj.x0 + obj.W; % initial noise support
            obj.uk = @(k) [0;0]; % fixed input function
            obj.T = 0.01;
            obj.G = 9.81;
            obj.L = 1;
            obj.m = 1;
            obj.J = obj.m * obj.L^2 / 12;
        end
        
        function x_ = f(obj,~,x,u)
			%F RHS of the ODE system
            [x1,x2,x3,x4] = unpack(x);
            z = zeros(size(x1));
            upart = [zeros(2);eye(2)]*u;
            if isa(x,'interval')
                upart = upart';
            end
            x_ = [x3;x4;z;z] + obj.n(x) + upart;
        end
        
        function n = n(obj,x)
            %N non-linear process function part
            %   x: current system state
            [x1,x2,x3,x4] = unpack(x);
            z = zeros(size(x1));
            s1 = sin(x1);
            s2 = sin(x2);
            s12 = sin(x1-x2);
            c12 = cos(x1-x2);
            n = 3./(16-9*c12.^2) .* (s12.*[z;z;-2*x4.^2-3*x3.^2.*c12;8*x3.^2+3*x4.^2.*c12]...
                +obj.G/obj.L.*[z;z;3*s2.*c12-6*s1;9*s1.*c12-8*s2]);
        end
        
        function z = h(obj,x)
            %H measurement function
            %   x: current system state
            [x1,x2] = unpack(x);
            z = [obj.L .* sin(x1);...
                -obj.L .* cos(x1);...
                obj.L .* sin(x1) + obj.L .* sin(x2);...
                -obj.L .* cos(x1) - obj.L .* cos(x2)];
%             z = [obj.L * sin(x1),...
%                 -obj.L * cos(x1),...
%                 obj.L * sin(x1) + obj.L * sin(x2),...
%                 -obj.L * cos(x1) - obj.L * cos(x2)];
%             if ~isa(x,'interval')
%                 z = z';
%             end
        end
        
        function o = ho(obj,z)
            %HO observator for measurment function
            %   z: measurement
            [z1,z2,z3,z4] = unpack(z);
            o = [atan((z4-z2)/(z3-z1));
                 atan(z2/z1)];
        end
        
        function E = energy(obj,x)
            %ENERGY calculate total energy of the system
            %   x: system states
            Ekin = obj.m * obj.L^2 / 6 * (4*x(3,:).^2 + x(4,:).^2 + 3*x(3,:).*x(4,:).*cos(x(1,:)-x(2,:)));
            Epot = obj.m * obj.G * obj.L / 2 * (4 - 3*cos(x(1,:)) - cos(x(2,:)));
            E = Ekin + Epot;
            E = E - E(1); % tare
        end
    end
end

