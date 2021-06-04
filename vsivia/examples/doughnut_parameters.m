%% Doughnut example
%
% This example is meant to show how to define a 2-dimensional inversion
% problem, in order to solve it using <../../html/vsivia.html vsivia>.
%
% We will consider the function _f_, defined as follows:
%
% $$ f : (x,y) \mapsto x^2 + y^2 + x \cdot y $$
% 
% The purpose will consist in determining an inner and an outer enclosures
% of the inverse image of _f_ on [1,2], using paving. The box [-3,3]x[-3,3]
% will be taken as an initial outer enclosure of the sought-after solution
% set.
%
% First of all, parameters for <../../html/vsivia.html vsivia> shall be
% specified in a class inherited from <../../html/vsivia_parameters.html
% vsivia_parameters>.

classdef doughnut_parameters < vsivia_parameters
    
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
        % _Y0_ gives the value to be inverted. Here, we consider [1,2].
        
        Y0 = [ 1 2 ] ;
        
        %%
        % _epsilon_ specifies the accuracy wanted for the result. A
        % relative _epsilon_ of 5 %, with respect to the size of _U0_, is
        % chosen: <../../html/vsivia.html vsivia> will not bisect any
        % boxes smaller than 5 % of the size of _U0_.
        
        epsilon = '1 %' ;
        
    end % properties
    
    
    %%
    % At last, the method _compute_ is overriden in order to define the
    % function _f_. Input variables of _compute_ appear here
    % separetely.
        
    methods (Static)
        
        function z = compute(x,y)
            
            z = x^2 + y^2 + x*y ;
            
        end
        
    end % methods
    
end % classdef
