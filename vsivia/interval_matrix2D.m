classdef interval_matrix2D < interval
    
    methods
    
        function obj = interval_matrix2D(varargin)
            
            obj = obj@interval(varargin{:}) ;
            
        end
        
        
        function r = mtimes(a,b)
            
           b = permute(b, [3 1 2]) ; 
           
           r = a .* b ;
           
           r = interval_matrix2D(permute(sum(r,2), [1 3 2])) ;
            
        end
        
        function r = ctranspose(a)  % r = r'        (transposes the interval array)
            r = interval_matrix2D(a.lower', a.upper') ;
        end % ctranspose
        
    end
    
end

