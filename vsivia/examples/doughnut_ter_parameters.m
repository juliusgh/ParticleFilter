%% Enhanced doughnut example
%
% Following the "doughnut example", this example shows how to define an
% inversion problem with two constraints. Let us consider thus:
% 
% $$ f : (x,y) \mapsto x^2 + y^2 $$
% $$ g : (x,y) \mapsto x + y $$
% 
% An interval _U_ will be considered as a solution of the problem if and
% only if it respects both constraints defined by _f_ and _g_.
%
% First of all, parameters for <../../html/vsivia.html vsivia> shall be
% specified in a class inherited from <../../html/vsivia_parameters.html
% vsivia_parameters>.

classdef doughnut_ter_parameters < vsivia_parameters
    
    properties
        
        %%
        % _algorithm_ tells <../../html/vsivia.html vsivia> the kind of
        % problem to be solved, here _inversion_.
        
        algorithm = 'inversion' ;
        
        %%
        % _U0_ indicates the initial box that is tested at first and
        % then possibly bisected. In this example, we consider the box
        % [-3,3]x[-3,3].
        
        U0 = [ -3 3 ; -3 3] ;
        
        %%
        % _Y0_ gives the value to be inverted. Here, we consider [0,4] with
        % respect to _f_ and [-10,1] with respect to _g_.
        
        Y0 = interval([ 0 2 ; -10 1 ])' ;
        
        %%
        % _epsilon_ specifies the accuracy wanted for the result. A
        % relative _epsilon_ of 1 %, with respect to the size of _U0_, is
        % chosen: <../../html/vsivia.html vsivia> will not bisect any
        % boxes smaller than 1 % of the size of _U0_.
        
        epsilon = '5%' ;
        
    end % properties
    
    
    %%
    % At last, the method _compute_ is overriden in order to define the
    % functions _f_ and _g_. Note that the result, namely _z_, returns an
    % _n_-by-2 array of intervals (_n_ being the length of vectors of
    % intervals _x_ and _y_), whose first column contains the images of
    % _x_ and _y_ through _f_, and whose second column contains the images
    % of _x_ and _y_ through _g_.
        
    methods (Static)
        
        function z = compute(x,y)
            
            z = [sqrt(x^2 + y^2,'p'), x + y] ;
            
        end
        
    end % methods
    
end % classdef
