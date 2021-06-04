classdef PS_vsivia < vsivia_parameters
    %PS_VSIVIA Set-Membership Filtering via Interval Analysis
    %   for 4-dim. ManipulatorSystem/PendulumSystem
    
    properties
        algorithm = 'inversion';
        U0;
        Y0;
        epsilon = '0.05%'; % accuracy
        maxBoxes = 5000;
        f; % measurement function handle
    end
    
    methods
        function obj = PS_vsivia(Y,U0)
            obj.U0 = U0;
            if isa(Y,'interval')
                obj.Y0 = Y;
            else
                obj.Y0 = interval(Y)';
            end
        end
        
        function z = compute(obj,x)
            z = obj.compute3(x);
%             if size(x,2) == 4
%                 z = obj.compute12(x);
%             else
%                 z = obj.compute1(x);
%             end
        end
        
        function z = compute12(obj,x)
            [x1,x2,x3,x4] = unpack(x);
            z = [sin(x1) + sin(x2);...
                -cos(x1) - cos(x2);...
                x3 .* cos(x1) + x4 .* cos(x2),...
                x3 .* sin(x1) + x4 .* sin(x2)];
        end
        
        function z = compute1(obj,x)
            [x1,x2] = unpack(x);
            z = [sin(x1) + sin(x2);...
                -cos(x1) - cos(x2)];
        end
        
        function z = compute3(obj,x)
            [x1,x2] = unpack(x);
            z = [sin(x1);...
                -cos(x1);...
                sin(x1) + sin(x2);...
                -cos(x1) - cos(x2)];
        end

    end % methods
end

