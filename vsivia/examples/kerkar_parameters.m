%% Modeling the intravenous injection of a drug into a subject
%
% This application presents a more advanced way to use
% <../../html/vsivia.html vsivia>, in comparison with
% <torus_parameters.html torus_parameters> and
% <doughnut_parameters.html doughnut_parameters>.
%
% The purpose of this application is to model, by a mathematical
% expression, the effect of a bolus intravenous injection of a drug into a
% subject. Concretely, we will aim at modeling the evolution of the blood
% concentration of the drug, using a function _f_ of the kind:
%
% $$ f : (a, \alpha, b, \beta, t) \mapsto a \cdot e^{- \alpha \cdot t} + b
% \cdot e^{- \beta \cdot t} $$
%
% where _t_ is the time elapsed since the injection. Enclosures of the four
% other parameters, that are constant parameters, will be computed using
% <../../html/vsivia.html vsivia>, from real measurements that give the
% value of _f_ for some known values of _t_. The measurements considered
% here are:
%
% $$ \begin{tabular}{|c||c|c|c|c|c|c|c|c|c|c|c|c|c|c|}
%    \hline
%    $t$ & 0.1 & 0.25 & 0.5 & 0.75 & 1 & 1.5 & 2 & 2.5 & 3 & 4 & 6 & 8 & 10 & 12 \\
%    \hline
%    $Y_0$ & 16.1 & 14.3 & 12 & 10.3 & 9 & 7.2 & 6.1 & 5.2 & 4.6 & 3.7 & 2.5 & 1.7 & 1.18 & 0.81 \\
%    \hline
%    \end{tabular} $$
%
% Let us set down $n$ the number of measurements, and _q_ the quadruplet of
% parameters we are looking for:
%
% $$ q = (a, \alpha, b, \beta) $$
%
% Let us introduce also _g_, defined as:
%
% $$ g : q \mapsto \left( \begin{tabular}{c} $f(q, t_1)$ \\ $f(q, t_2)$ \\ ... \\ $f(q, t_n)$ \end{tabular} \right) $$
%
% In this way, the problem of determining _q_ can be seen as the inversion
% of _Y0_ through _g_. Hence, enclosures of _q_ may be found using
% <../../html/vsivia.html vsivia>.
%
% The following class shows how to configure
% <../../html/vsivia.html vsivia> in order to solve this problem.

%%
% To begin with, a new class inheriting from
% <../../html/vsivia_parameters.html vsivia_parameters> shall be created.

classdef kerkar_parameters < vsivia_parameters
    
    properties
        
        %%
        % Its _algorithm_ property shall be _inversion_.
        
        algorithm = 'inversion' ;
        
        %%
        % An initial guess about a box containing _q_ shall be specified.
        
        U0 = [8 20 ; -6 6] ;
        
        %%
        % Then the image to be inverted, namely the measurements.
        
       Y0 = interval([11.3623,11.9576; 6.4199,6.9998; 4.6245,5.2251; 5.5761,6.1794; 9.7981,10.3528])';
       
        %%
        % Then the accuracy parameter, _epsilon_, that is set to 10 %, with
        % respect to the size of _U0_.
        
        epsilon = '5 %' ;
        
    end % properties
    
    %%
    % In addition to the "core parameters" introduced by
    % <../../html/vsivia_parameters.html vsivia_parameters>, the time vector
    % is introduced as an additionnal property.
    
    properties (Constant)
        
        t = [0 20 40 60 80] ;
        
    end % constant properties
    
    %%
    % Finally, the function to be inverted (_i.e._ _g_) is defined. In
    % comparison with <doughnut_parameters.html doughnut_parameters>
    % and <torus_parameters.html torus_parameters>, the input
    % parameter, namely _q_, is processed as a single parameter, and not as
    % four separate parameters.

    methods (Static)
        
        function y = compute(p)
            
            t2 = kerkar_parameters.t ;
            
            y = p(:,1)*cos(t2)+sqrt((p(:,2))^2 -(p(:,1)*sin(t2))^2);
            
        end
        
    end % methods
    
end
