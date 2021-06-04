classdef ManipulatorSystemDirect < ManipulatorSystem
    %MANIPULATORSYSTEM "Ebener Manipulator"
    %   Simple mechanical test system with direct state measurement
    
    methods
        function obj = ManipulatorSystemDirect()
            %MANIPULATORSYSTEMDIRECT Construct an instance of this class
            %   Simple mechanical test system with direct state measurement
            obj = obj@ManipulatorSystem();
            obj.W = HyperRect([-0.3 0.3],[-0.3 0.3],[-0.3 0.3],[-0.3 0.3]); % measurement noise support
            obj.Wsim = obj.W;
            obj.x0 = [pi;pi;0;0]; % initial condition
            obj.X0 = obj.x0 + obj.W; % initial noise support
        end
        
        function z = h(~,x)
            %H measurement function
            %   x: current system state
            %z = [eye(2) zeros(2)]*x;
            z = x;
        end
    end
end

