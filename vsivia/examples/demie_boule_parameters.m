%% Hollow Hemisphere example
%
% In this example, we will aim at defining a 3-dimensional inversion
% problem, in order to solve it using <../../html/vsivia.html vsivia>.
%
% The following function _f_ will be considered:
%
% $$ f : (x,y,z) \mapsto ((x^2 + y^2+z^2), x) $$
%
% together with the constraints $ Y0[1,1) < f_1(x,y,z) < Y0[1,2) $
% together with the constraints $ Y0[2,1) < f_2(x,y,z) < Y0[2,2) $
% 
% <../../html/vsivia.html vsivia> is used to determine an inner and an
% outer enclosures of the inverse image of _f_ on [Y0[:,1),Y0[:, 2)], using paving. 
% From a graphical point of view, we are looking for the intersection of 
%                  - the semi-space $x>0$,
%                  - the space between two spheres whose radii are 
%                  resp.  $ \sqrt(Y0[1,1))$ and   $ \sqrt(Y0[1,2))$. 
% We set the initial box to
% [-8,8]x[-8,8]x[-8,8].
%
% At first, we shall create a class inheriting from
% <../../html/vsivia_parameters.html vsivia_parameters>

classdef demie_boule_parameters < vsivia_parameters
    
    properties
        
        %%
        % We indicate the kind of problem: _inversion_.
        
        algorithm = 'inversion' ;

        %%
        % Then the initial box.
        
        U0 = [-8 8 ; -8 8 ; -8 8] ;
        
        %%
        % Then the interval image to be inverted.
        
        Y0 = interval([1 4 ; 0 10 ])' ;
        
        %%
        % Then the accuracy parameter. We consider here an absolute epsilon of 0.5,
        % which means that boxes whose larger component is smaller than 0.5
        % will not be bisected.
        
        epsilon = .1 ;
        
        
    end % properties
    
    %%
    % At last, the function _f_ is defined.
    
    methods (Static)
        
        function z = compute(x,y,z)
                        
            z = [x^2+y^2+z^2, x] ;
            
        end
        
    end % methods
    
end

