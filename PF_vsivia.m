classdef PF_vsivia < vsivia_parameters
    %PF_VSIVIA Parameters for vsivia
    
    properties
        algorithm = 'inversion';
        U0;
        Y0;
        epsilon = '2%'; % accuracy
        maxBoxes = 30000; % maximum number of boxes or NaN
        S;
        f;
        intersect;
    end
    
    methods
        function obj = PF_vsivia(f,Y,U0,S,intersect)
            obj.U0 = U0;
            obj.f = f;
            obj.S = S;
            if ~isa(Y,'interval')
                Y = interval(Y)';
            end
            obj.intersect = intersect;
            if intersect
                obj.Y0 = [interval(Y);interval([1 1])];
            else
                obj.Y0 = interval([1 1]);
            end
        end
        
        function z = compute(obj,x)
            if obj.intersect
                z1 = obj.f(x);
                in = max(isinV(x,obj.S),[],2);
                isects = ~min(isoutV(x,obj.S),[],2);
                z = [z1;interval([in isects])];
            else
                %z = z1;
                eps = 0.001; % eps for initial mincing
                in = min(bsxfun(@lt, width(x), eps), [], 2);
                o = ones(size(in));
                z = interval(in,o);
            end
        end
        
        function z = f0(obj,x,eps)
            z = min(bsxfun(@lt, width(x), eps), [], 2);
        end


    end % methods
end

