classdef Grid
    %GRID Simple grid implementation
    %   Represents a hyperrectangle grid with constant edge-length
    
    properties
        eps % edge-length
        lb % lower bound of grid for each axis
        ub % upper bound of grid for each axis
        counts % cell counts for each axis
        dim % grid dimensions
        cells % cell array
    end
    
    methods
        function obj = Grid(X,eps)
            %GRID Construct an instance of this class
            %   Create grid with edge-length eps around X
            obj.eps = eps;
            obj.dim = size(X,1);
            obj.lb = NaN(1,obj.dim);
            obj.ub = NaN(1,obj.dim);
            obj.counts = NaN(1,obj.dim);
            for i = 1:obj.dim % loop through every dimension of X
                obj.lb(i) = min(X(i,:)); % determine lower bound
                obj.ub(i) = max(X(i,:)); % determine upper bound
                obj.counts(i) = floor((obj.ub(i) - obj.lb(i)) / obj.eps) + 1;
            end
            c = num2cell(obj.counts);
            obj.cells = zeros(c{:});
        end
        
        function free = freeCell(obj,x)
            %FREECELL Check if cell of x is free
            %   Returns logical 1 or logical 0
            if isnan(x)
                free = false;
            else
                c = obj.getCoords(x);
                value = obj.cells(c{:});
                free = ~value;
            end
        end
        
        function obj = occupyCell(obj,x)
            %OCCUPYCELL Mark cell of x as occupied
            %   Returns new grid
            if ~isnan(x)
                c = obj.getCoords(x);
                obj.cells(c{:}) = 1;
            end
        end
        
        function coords = getCoords(obj,x)
            %GETCELL Get coordinates of the cell around x
            %   Returns coordinates of cell
            c = NaN(obj.dim,1);
            for i = 1:obj.dim
                c(i) = floor((x(i) - obj.lb(i)) / obj.eps) + 1;
            end
            coords = num2cell(c);
        end
    end
end

