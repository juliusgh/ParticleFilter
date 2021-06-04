classdef HyperRect
    %HYPERRECT Simple implementation of a hyperrectangle set
    %   Provides random sample and element check
    
    properties
        bounds % bounds of set
        dim % dimensions of set
    end
    
    methods
        function obj = HyperRect(varargin)
            %HYPERRECT Construct an instance of this class
            %   Input: [lb1 ub1],...,[lbn lbn]
            if nargin > 0
                obj.dim = nargin;
                if nargin == 1
                    obj.dim = size(varargin{1},1);
                    varargin = mat2cell(varargin{1},ones(1,obj.dim),2);
                end
                obj.bounds = NaN(obj.dim,2);
                for i = 1:obj.dim
                    obj.bounds(i,:) = varargin{i};
                end
            end
        end
        
        function x = sample(obj,N)
            %SAMPLE Generates N random samples from set
            %   Returns N random samples (default: N=1)
            if nargin < 2
                N = 1;
            end
            lower = obj.bounds(:,1);
            upper = obj.bounds(:,2);
            x = lower + (upper - lower) .* rand(obj.dim,N);
        end
        
        function x = sampleBorder(obj,N)
            %SAMPLE Generates N random samples from border of set
            %   Returns N random samples (default: N=1)
            if nargin < 2
                N = 1;
            end
            B = obj.bounds';
            ind = randi([1 2],obj.dim,N)' + (0:2:2*(obj.dim-1));
            x = B(ind)';
        end
        
        function x = sampleCC(obj,N)
            %SAMPLE Generates N random samples with counter-comonotonicity
            %   Returns N random samples (default: N=1)
            if nargin < 2
                N = 1;
            end
            lower = obj.bounds(:,1);
            upper = obj.bounds(:,2);
            r = (-1).^(0:obj.dim-1) * (1 - 2*unifrnd(0,1).^2);
            x = lower + (upper - lower) .* r';
            x = 0.03 * [1; -1] * ( 1 - 2*unifrnd(0,1).^2 );
        end
        
        function r = toInterval(obj)
            %TOINTERVAL Return equivalent `interval` object
            %   `interval` is a class from the VSIVIA framework
            r = interval(obj.bounds(:,1)',obj.bounds(:,2)');
        end
        
        function w = width(obj)
            %WIDTH Return the width of the set
            %   width is defined as the largest edge-length
            lower = obj.bounds(:,1);
            upper = obj.bounds(:,2);
            w = max(abs(upper - lower));
        end
        
        function pa = principalAxis(obj)
            %PRINCIPALAXIS Return the principal axis of the set
            %   principal axis is defined as the axis with the largest edge-length
            lower = obj.bounds(:,1);
            upper = obj.bounds(:,2);
            [~,pa] = max(upper - lower);
        end
        
        function in = contains(obj,x)
            %CONTAINS Check if set contains element x
            %   Returns logical 1 or logical 0
            lower = obj.bounds(:,1);
            upper = obj.bounds(:,2);
            if isnumeric(x)
                in = min(x >= lower & x <= upper);
            else
                xlower = x.bounds(:,1);
                xupper = x.bounds(:,2);
                in = min(xlower >= lower & xupper <= upper);
            end                
        end
        
        function in = intersects(obj,x)
            %CONTAINS Check if set contains element or set x
            %   Returns logical 1 or logical 0
            lower = obj.bounds(:,1);
            upper = obj.bounds(:,2);
            if isnumeric(x)
                in = obj.contains(x);
            else
                xlower = x.bounds(:,1);
                xupper = x.bounds(:,2);
                in = min((xlower >= lower & xlower <= upper) | (xupper >= lower & xupper <= upper)) | ...
                     min((lower >= xlower & lower <= xupper) | (upper >= xlower & upper <= xupper));
            end                
        end
        
        function [in,out] = inclusionTest(obj,x)
            xb = [x.bounds];
            xlb = xb(:,1:2:end);
            xub = xb(:,2:2:end);
            b = repmat(obj.bounds,1,numel(x));
            lb = b(:,1:2:end);
            ub = b(:,2:2:end);
            in = min((xlb >= lb)&(xub <= ub),[],1);
            out = ~(min((xlb >= lb & xlb <= ub) | (xub >= lb & xub <= ub),[],1) | ...
                    min((lb >= xlb & lb <= xub) | (ub >= xlb & ub <= xub),[],1));
        end
        
        function paving = SIVIA(obj,f,x,eps)
        %SIVIA Construct a paving of the inverse set given the function f
        %   x: prior set, eps: accuracy, maxk: maximum iterations
            if nargin < 4
                eps = 0.02;
            end
            % zweites f bzw. f erweitern
            % f1 = [1 1]*x
            % f2 = foreach u: u und x disjunkt
            % f2 = ~min(u und x disjunkt) über alle u
            % TODO: schnitt mit and überladen
            % u: Boxen aus letztem Zeitschritt propagiert
            % obj erweitern um 2. Dimension: [1,1]
            % alternativ: upper+lower bound = 1,
            stack = x;
            Kin = HyperRect.empty();
            Ku = HyperRect.empty();
            time1 = 0;
            time2 = 0;
            while ~isempty(stack)
                t1 = tic;
                y = arrayfun(f,stack);
                time1 = time1 + toc(t1);
                t2 = tic;
                [in,out] = obj.inclusionTest(y);
                u = (~in)&(~out);
                Kin = [Kin stack(in)];
                U = stack(u);
                ue = arrayfun(@(x) x.width() < eps,U);
                Ku = [Ku U(ue)];
                b = arrayfun(@(x) x.bisect(),U(~ue),'UniformOutput',false);
                stack = [b{:}];
                time2 = time2 + toc(t2);
            end
            disp(time1);
            disp(time2);
            paving = Paving(Kin,Ku);
        end
        
        function paving = classicSIVIA(obj,f,x,eps,maxk)
        %SIVIA Construct a paving of the inverse set given the function f
        %   x: prior set, eps: accuracy, maxk: maximum iterations
            if nargin < 5
                maxk = 1e5;
                if nargin < 4
                    eps = 0.02;
                end
            end
            % zweites f bzw. f erweitern
            % f1 = [1 1]*x
            % f2 = foreach u: u und x disjunkt
            % f2 = ~min(u und x disjunkt) über alle u
            % TODO: schnitt mit and überladen
            % u: Boxen aus letztem Zeitschritt propagiert
            % obj erweitern um 2. Dimension: [1,1]
            % alternativ: upper+lower bound = 1,
            stack = [];
            Kin = HyperRect.empty();
            Ki = HyperRect.empty();
            for k = 1:maxk
                y = f(x);
                if obj.contains(y)
                     Kin(end+1) = x;
                elseif ~obj.intersects(y)
                    % do nothing...
                elseif x.width() <= eps
                    Ki(end+1) = x;
                else
                    stack = [stack x.bisect()];
                end
                if numel(stack) > 0
                    x = stack(1);
                    stack = stack(2:end);
                else
                    break;
                end
            end
            if k == maxk
                disp('reached maximum iterations!');
            end
            paving = Paving(Kin,Ki);
        end
        
        function parts = divide(obj,divisions)
            % split HyperRect along every axis in `divisions` parts
            n = divisions^(obj.dim);
            parts = HyperRect.empty(0,n);
            points = NaN(obj.dim,divisions);
            steps = NaN(obj.dim,1);
            lower = obj.bounds(:,1);
            upper = obj.bounds(:,2);
            for d = 1:obj.dim % vectorization?
                ls = linspace(lower(d),upper(d),divisions+1);
                steps(d) = (upper(d)-lower(d))/divisions;
                points(d,:) = ls(1:end-1);
            end
            % ndgrid over all dimensions
            grids = mat2cell(points,ones(1,obj.dim));
            [grid{1:obj.dim}] = ndgrid(grids{:});
            % construct parts
            for i = 1:n
                params = cell(1,obj.dim);
                for d = 1:obj.dim
                    lb = grid{d};
                    ub = lb + steps(d);
                    params{d} = [lb(i) ub(i)];
                end
                parts(i) = HyperRect(params{:});
            end
        end
        
        function parts = bisect(obj)
            % split HyperRect along principal axis in 2 parts
            a = obj;
            b = obj;
            pa = obj.principalAxis();
            a.bounds(pa,2) = sum(a.bounds(pa,:)) / 2;
            b.bounds(pa,1) = sum(b.bounds(pa,:)) / 2;
            parts = [a b];
        end
        
        function v = vertices(obj)
            if obj.dim == 2
                v = [obj.bounds(1,1) obj.bounds(2,1); obj.bounds(1,1) obj.bounds(2,2);
                     obj.bounds(1,2) obj.bounds(2,1); obj.bounds(1,2) obj.bounds(2,2)];
            else
                v = [];
                disp('not supported yet');
            end
        end
        
        function p = polygon(obj)
            if obj.dim == 2
                p = polyshape([obj.bounds(1,1) obj.bounds(1,1) obj.bounds(1,2) obj.bounds(1,2)],...
                              [obj.bounds(2,2) obj.bounds(2,1) obj.bounds(2,1) obj.bounds(2,2)]);
            else
                p = polyshape();
                disp('not supported yet');
            end
        end
        
        function plot(obj, color)
            if nargin < 2
                color = 'g';
            end
            if obj.dim == 2
                rectangle('Position',[obj.bounds(1) obj.bounds(2) obj.bounds(3)-obj.bounds(1) obj.bounds(4)-obj.bounds(2)],...
                    'EdgeColor', color);
                %axis([0 5 0 5]);
            else
                disp('not supported yet');
            end
        end
        
        function r = get(obj,i)
            r = HyperRect(obj.bounds(i,:));
        end
        
        function obj = set(obj,i,r)
            obj.bounds(i,:) = r.bounds(1,:);
        end
        
%         function varargout = subsref(obj,s)
%             [varargout{1:nargout}] = builtin('subsref',obj,s);
%         end
        
        function r = vertcat(varargin)
            c = cellfun(@(x) x.bounds,varargin,'UniformOutput',false);
            r = HyperRect(vertcat(c{:}));
        end
        
        function r = plus(r1,r2)
            if isnumeric(r1) && isnumeric(r2)
                r = r1 + r2;
            elseif isnumeric(r1)
                r = HyperRect(r1 + r2.bounds);
            elseif isnumeric(r2)
                r = HyperRect(r1.bounds + r2);
            else
                r = HyperRect(r1.bounds + r2.bounds);
            end
        end
        
        function r = minus(r1,r2)
            if isnumeric(r1) && isnumeric(r2)
                r = r1 - r2;
            elseif isnumeric(r1)
                r = HyperRect(r1 - r2.bounds);
            elseif isnumeric(r2)
                r = HyperRect(r1.bounds - r2);
            else
                r = HyperRect(r1.bounds - r2.bounds);
            end
        end
        
        function r = mtimes(r1,r2)
            if isnumeric(r1) || isnumeric(r2)
                if isnumeric(r1)
                    a = r2;
                    b = r1;
                else
                    a = r1;
                    b = r2;
                end
                [n,m] = size(b);
                c = cell(n,1);
                if m == a.dim
                    for i = 1:n
                        c{i} = HyperRect([0 0]);
                        for j = 1:m
                            c{i} = c{i} + b(i,j) .* a.get(j);
                        end
                    end
                    r = vertcat(c{:});
                else
                    r = b .* a;
                end
            else
                r = r1 .* r2;
            end
        end
        
        function r = times(r1,r2)
            if isnumeric(r1) && isnumeric(r2)
                r = r1 * r2;
            elseif isnumeric(r1)
                r = HyperRect(r1 * r2.bounds);
            elseif isnumeric(r2)
                r = HyperRect(r1.bounds * r2);
            else
                c = [r1.bounds(1)*r2.bounds(1),r1.bounds(1)*r2.bounds(2),r1.bounds(2)*r2.bounds(1),r1.bounds(2)*r2.bounds(2)];
                r = HyperRect([min(c), max(c)]);
            end
        end
        
        function r = mpower(a,b)
            r = a .^ b;
        end
        
        function r = power(a,b)
           if b == 1
               r = a;
           elseif b == 2
               r = a .* a;
           else
               disp('not implemented yet');
           end
        end
    end
end

