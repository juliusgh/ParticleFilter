%% Solving an IVGTT problem using SIVIA
%
% Minimal models of plasma glucose and insulin kinetics have been developed
% and used by Bergman and co-workers since the 1970's to investigate
% glucose metabolism _in vivo_ in physiological, pathological and
% epidemiological studies from a frequently-sampled intravenous glucose
% tolerance (FSIGT) test, _i.e._ standard intravenous glucose tolerance test
% (IVGTT). Amongst others indices, the glucose and insulin minimal models
% allow us to characterize the FSIGT test data in terms of four metabolic
% indices:
%
% * _Si_ : insulin sensitivity, the dependence of fractional glucose
% disappearance on plasma insulin
%
% * _Sg_ : glucose effectiveness, the ability of glucose _per se_ to
% suppress endogenous glucose production and stimulate glucose uptake
%
% * _p1_ : first phase pancreatic responsivity, a measurement of the size
% of the first peak in plasma insulin due to the glucose injection
%
% * _p2_ : second phase pancreatic responsivity, a measurement of the size
% of the second peak of plasma insulin that follows the first peak and the
% refractory period


%%
% In a typical FSIGT test, blood samples are taken from a fasting subject
% at regular intervals of time, following a single intravenous injection of
% glucose (_t_ = 0). The blood samples are then analysed for glucose and
% insulin contents. The table below gives a typical response from a healthy
% subject to a single intravenous injection of glucose:
%
% $$ \begin{tabular}{|l||c|c|c|c|c|c|c|c|c|c|c|c|}
%    \hline
%    $t$ (min) & 0 & 2 & 4 & 6 & 8 & 10 & 12 & 14 & 16 & 19 & 22 & 27 \\
%    \hline
%    $G$ (mg/dL) & 92 & 350 & 287 & 251 & 240 & 216 & 211 & 205 & 196 & 192 & 172 & 163 \\
%    \hline
%    $I$ (mU/mL) & 11 & 26 & 130 & 85 & 51 & 49 & 45 & 41 & 35 & 30 & 30 & 27 \\
%    \hline
%    \end{tabular} $$
%
% $$ \begin{tabular}{|l||c|c|c|c|c|c|c|c|c|c|c|c|}
%    \hline
%    $t$ (min) & 32 & 42 & 52 & 62 & 72 & 82 & 92 & 102 & 122 & 142 & 162 & 182 \\
%    \hline
%    $G$ (mg/dL) & 142 & 124 & 105 & 92 & 84 & 77 & 82 & 81 & 82 & 82 & 85 & 90 \\
%    \hline
%    $I$ (mU/mL) & 30 & 22 & 15 & 15 & 11 & 10 & 8 & 11 & 7 & 8 & 8 & 7 \\
%    \hline
%    \end{tabular} $$

%%
% Denoting _Gb_ and _Ib_ respectively the glucose and insulin steady-state
% concentrations (also referred to as basal concentrations, usually
% determined by averaging the last IVGTT measurements), Bergman minimal
% model describes the evolution of the blood glucose concentration by the
% differential system hereafter:
%
% $$ \left\lbrace \begin{array}{ccl}
%    \dot{X} & = & p_2 \cdot (S_i \cdot (I - I_b) - X) \\
%    \dot{G} & = & - (S_g + X) \cdot G + S_g \cdot G_b \\
%    X(0) & = & 0 \\
%    G(0) & = & G_0
%    \end{array} \right. $$

%%
% Hence, the purpose of this application is to find the possible values of
% the parameters introduced in Bergman minimal model (_Si_, _Sg_, _p2_ and
% _G0_), in order to make the model fit the measurements above.

%%
% The approach is quite similar to that followed in
% <drugs_parameters.html drugs_parameters> : the model is
% computed for each _t_ of the measurement table, and its results compared
% to the corresponding measurements. A quadruplet (_Si_, _Sg_, _p2_, _G0_)
% is considered as a solution of the problem if and only if it makes the
% model fit every measurement.

%%
% Using interval arithmetics, for this application, makes it possible to
% find every quadruplet (_Si_, _Sg_, _p2_, _G0_), such as Bergman minimal
% model fits a set of data subject to measurements errors. Here, a 2 %
% relative error on _G_, as well as a relative 3 % error on _I_, are
% considered. As done in <drugs_parameters.html drugs_parameters>,
% determining the possible values of (_Si_, _Sg_, _p2_, _G0_) can be seen
% as a set inversion problem, and may be performed using
% <../../html/vsivia.html vsivia>.

%%
% However, in comparison with <drugs_parameters.html drugs_parameters>,
% the considered function is not trivial to compute, since no explicit
% expression of it is known. It has thus to be integrated. Here, we will
% use a simple explicit Euler method to achieve this. For a given time
% interval _dt_ :
%
% $$ \left\lbrace \begin{array}{ccl}
%    X(t + dt) & = & X(t) + dt \cdot p_2 \cdot (S_i \cdot (I(t) - I_b) - X(t)) \\
%    G(t + dt) & = & G(t) + dt \cdot (-(S_g + X(t)) \cdot G(t) + S_g \cdot G_b)
%    \end{array} \right. $$

%%
% Computing _X_ and _G_ with this method requires to have measurements of
% _I_ at any time (at least every _dt_), which is unlikely to occur. Even
% though it introduces an estimation error, interpolation is thus used to
% solve this issue.
%
% An other source of estimation error is caused by overbounding phenomena,
% that may happen in interval computations when a same interval occurs
% several times in a same expression. Here, this problem is raised for both
% _X_ and _G_.

%%
% Let us change the basis of the problem by introducing _p3_, such as :
%
% $$ p_3 = p_2 \cdot S_i $$
%
% The system may be rewritten as follows:
%
% $$ \left\lbrace \begin{array}{ccl}
%    X(t + dt) & = & X(t) \cdot (1 - dt \cdot p_2) + dt \cdot p_3 \cdot (I(t) - I_b) \\
%    G(t + dt) & = & G(t) \cdot (1 - dt \cdot (S_g + X(t))) + dt \cdot S_g \cdot G_b
%    \end{array} \right. $$

%%
% The only remaining issue concerns the bi-occurence of _Sg_ in the
% expression of _G_. It can be solved using the *-partially optimal
% coercion theorem, according to which the modality of _Sg_ can be switched
% in the following way:
%
% $$ \left\lbrace \begin{tabular}{ccccl}
%     $G(t + dt)$ & $=$ & $G(t) \cdot (1 - dt \cdot (S_G + X(t))) + dt \cdot S_G^\ast \cdot G_b$ && if $G \geq G_b$ \\
%     $G(t + dt)$ & $=$ & $G(t) \cdot (1 - dt \cdot (S_G^\ast + X(t))) + dt \cdot S_G \cdot G_b$ && if $G \leq G_b$ \\
%     $G(t + dt)$ & $=$ & $G(t) \cdot (1 - dt \cdot (S_G + X(t))) + dt \cdot S_G \cdot G_b$ && otherwise\\
%    \end{tabular} \right. $$

%%
% To conclude, the IVGTT problem can be solved as an inversion problem. In
% order to reduce overbounding phenomena in the computations, modal
% interval arithmetics is used, and an intermediate variable, _p3_, has
% been introduced to replace _Si_ in the computations. Once a solution of
% the problem has been found, _Si_ can be simply determined from _p2_ and
% _p3_.

%%
% In comparison with <doughnut_parameters.html doughnut_parameters>,
% <torus_parameters.html torus_parameters>, and <drugs_parameters.html
% drugs_parameters>, the class used to define the parameters of this IVGTT
% application for <../../html/vsivia.html vsivia> involves additionnal
% processing, in particular importing measurements from a file,
% interpolating the measurements of _I_, and performing the integration of
% _X_ and _G_. However, the way to write the class remains the same, since
% the class shall still derivate from <../../html/vsivia_parameters.html
% vsivia_parameters>.

classdef IVGTT_parameters < vsivia_parameters
    
    properties
        
        %%
        % The problem is an _inversion_ problem.
        
        algorithm = 'inversion' ;
        
        %%
        % The following box is taken as an initial guess of an enclosure of
        % the solutions (_p1_, _p2_, _p3_, _G0_).
        
        U0 = [.01 .04 ; .01 .04 ; .5e-5 5e-5 ; 200 350] ;
        
        %%
        % Declaration of the measurement vector, _Y0_, that is initialized
        % in the constructor.
        
        Y0 ;
        
        %%
        % The accuracy parameter, upon which the length of the computation
        % depends.
        
        epsilon = '2^-4 rel';
        
        %%
        % A tolerance parameter telling <../../html/vsivia.html vsivia>
        % that up to 6 mismatched measurements are acceptable for a
        % solution of the problem. This tolerance is considered because the
        % model cannot keep up with the significant variations of the
        % earliest measurements.
        
        tol_out = 6 ;
        
        %%
        % The time vector for the measurements of _G_ (loaded later).
        
        tG ;
        
        %%
        % The time vector for the measurements of _I_ (loaded later).
        
        tI ;
        
        %%
        % The measurements of _G_, taking the measurement uncertainties
        % into account (computed later).
        
        Gerr ;
        
        %%
        % The measurements of _I_, taking the measurement uncertainties
        % into account (computed later).
        
        Ierr ;
        
    end % properties
    
    
    methods
        
        %%
        % The measurements of _G_ and _I_ are loaded when instantiating the
        % class, _i.e._ in the constructor hereafter. The measurement
        % errors are also computed in this method. The measurement vector
        % for <../../html/vsivia.html vsivia> corresponds to the
        % measurements of _G_, taking the measurement errors into account.
        
        function obj = IVGTT_parameters
            
            load IVGTT_data ;
            
            obj.tG = t' ;
            obj.tI = t' ;
            
            G0 = interval(G, [], 0)' ;
            I0 = interval(I, [], 0)' ;
            
            err = max(2.5, .02 * G') ;
            
            obj.Gerr = G0 + interval(-err, err) ;
            obj.Ierr = [ .97 1.03] * I0 ;
            
            obj.Y0 = obj.Gerr ;
            
        end % IVGTT_parameters
        
        
        %%
        % Lastly, the inversion function is defined. Note here that it is
        % declared as an instance method and not as a static method, in
        % order to use the values loaded and computed in the class
        % constructor, that are stored within instance properties.
        %
        %
        
        function y = compute(this, p1, p2, p3, g0)
            
            % Loading of the time and measurement vectors...
            
            tG_ = this.tG ;
            tI_ = this.tI ;
            
            G = this.Gerr ;
            I = this.Ierr ;
            
            % Basal glucose and insulin are taken as the final measurements
            % of G and I (when they become almost constant).
            
            Gb = G(end) ;
            Ib = I(end) ;
            
            % Time step and time vector for the numerical integration.
            
            dt = 1 ;
            
            it = 0:dt:tI_(end) ;
            
            % Interpolation of the measurements of I...
           
            Ii = interval(interp1(tI_,I.lower,it,'cubic'), interp1(tI_,I.upper,it,'cubic')) ;
            
            % Numerical integration to compute X and G...
            
            n = size(p1,1) ;
            
            nit = numel(it) ;
            
            Xe = interval(zeros(n,nit), [], 0) ;
            Ge = interval(zeros(n,nit), [], 0) ;
            
            Ge(:,1) = g0 ;
            
            for k=1:nit-1
                
                Xe(:,k+1) = Xe(:,k) * (1 - dt*p2) + dt*p3*(Ii(k)-Ib) ;
                
                geGgb = Ge(:,k) >= Gb ;
                geLgb = Ge(:,k) <= Gb ;
                geIgb = ~(geGgb|geLgb) ;
                
                Ge(geGgb,k+1) = Ge(geGgb,k) * (1 - dt*(Xe(geGgb,k) + p1(geGgb))) + dt * p1(geGgb).' * Gb ;
                
                Ge(geLgb,k+1) = Ge(geLgb,k) * (1 - dt*(Xe(geLgb,k) + p1(geLgb).')) + dt * p1(geLgb) * Gb ;
                
                Ge(geIgb,k+1) = Ge(geIgb,k) * (1 - dt*(Xe(geIgb,k) + p1(geIgb))) + dt * p1(geIgb) * Gb ;
                
            end
            
            % Interpolation to obtain the results of the model when the
            % measurements of G have been performed.
            
            y = interval(interp1(it,Ge.lower',tG_,'cubic'), interp1(it,Ge.upper',tG_,'cubic')) ;
            
            if n > 1
                y = y' ;    % Dimension 1: boxes ; dimension 2: interpolated data
            end
            
        end % compute
        
    end % methods
    
end

