%% Using VSIVIA
%
%%
% <vsivia.html vsivia> is part of a toolbox that comprises several
% components: <vsivia.html vsivia> itself, a vectorial interval library,
% a contractor library, as well as some other tool functions, in particular
% <draw_boxes.html draw_boxes>, that is meant to plot sets of boxes on a
% figure. The pieces of code hereafter, based on the examples provided with
% <vsivia.html vsivia>, show how to implement these features.

clear all ;
close all ;

%% Solving a 2-D inversion problem
% 
% <../examples/html/doughnut_parameters.html doughnut_parameters> defines
% the following problem: inverting the function _f_ on [1,2], _f_ being:
%
% $$ f : (x,y) \mapsto x^2 + y^2 + x \cdot y $$
%
% <draw_boxes.html draw_boxes> is used in order to plot the results of
% <vsivia.html vsivia>. The three first parameters specify the three
% vectors of boxes to be plotted, which are displayed respectively in red,
% yellow and blue. The last parameters indicate the dimensions associated
% to the X, Y and Z-axes. Their number depends upon the number of
% dimensions of the considered problem.
%
% On the figure below, red boxes are sets whose image through _f_ lies
% between 1 and 2: the set of these boxes is an inner enclosure of the
% solution set of the problem. On the opposite, blue boxes have been proven
% to contain no solution at all. Between those red and blue sets, yellow
% boxes are boxes whose images through _f_ are not included in [1,2], but
% that have an non empty intersection with in [1,2]. Thus, the set of red
% and yellow boxes is an outer enclosure of the solution set.

%% doughnut 1
% [S,E,N] = vsivia_lite(doughnut_parameters) ;   
% draw_boxes(S,E,N,1,2) ;


% %%doughnut 2
[S,E,N] = vsivia(doughnut_ter_parameters) ;
draw_boxes(S,E,N,1,2) ;


%% Solving a 3-D inversion problem
%
% <../examples/html/torus_parameters.html torus_parameters> defines a
% similar problem: inverting _g_ on [1,4], where:
%
% $$ g : (x,y,z) \mapsto (R - \sqrt{x^2 + y^2})^2 + z^2 $$
%
% In this example, _R_ is set to 5. As in the previous application, the
% results are plotted using <draw_boxes.html draw_boxes>. When called with
% 6 parameters, the figure is in 3 dimensions and, for a good readability,
% only the centers of the boxes contained in _S_ and _E_ are drawn using
% red and yellow spheres.

% [S,E,N] = vsivia(torus_parameters) ;
% draw_boxes(S,E,N,1,2,3) ;


%% Mechanical problems
% C is the center of the wheel
% the matrix F defines the normal vectors 
%            F(:,i) is the ith normal
% Nf is the number of planes that define the cone ()
%   
% 
% % Problem 1
% [S,E,N] = vsivia(demie_boule_parameters) ;
%     
% draw_boxes(S,E,N,1,2,3) ;
% axis equal;

% % % Problem 2
% crsa = contact_roue_sol_a;
% [S,E,N] = vsivia(crsa) ;
% draw_boxes(S,E,N,1,2,3) ;
% axis equal;


%% Solving a 4-D data-fitting problem
%
% In this third example, a data-fitting problem is solved as an inversion
% problem. As defined in <../examples/html/drugs_parameters.html drugs_parameters>,
% an _n_-element time vector _t_ as well as an _n_-element measurement
% vector _Y0_ are considered, the purpose being to find the boxes _x_ such as:
% 
% $$ f(t_i, x) \; \in \; Y_0^i ~~ (\forall i \in [1,n]) $$
%
% _f_ being:
%
% $$ f : (x_1, x_2, x_3, x_4, t) \mapsto x_1 \cdot e^{-x_2 \cdot t}
%           + x_3 \cdot e^{-x_4 \cdot t} $$
%
% On the figure below, measurements have been drawn in black, using error
% bars in order plot measurement uncertainties, _i.e._ the lower and upper
% possible values for _Y0_. From the parameters defined by
% <../examples/html/drugs_parameters.html drugs_parameters>, <vsivia.html
% vsivia> determines an inner and an outer enclosures of the solution set,
% that corresponds to the situation in which _f_ matches every measurement.

%drugs = kerkar_parameters() ;

% drugs = drugs_parameters() ;
% 
% [S,E,N] = vsivia_lite(drugs) ;
% 
% t = drugs.t ;
% 
% Y0 = drugs.Y0 ;
% 
% figure ; 
% 
% hold on ;
% 
% errorbar(t, .5*(Y0.lower+Y0.upper), .5*(Y0.upper-Y0.lower), 'kx')
% 
% if ~isempty(S)
%     YS = join(drugs.compute(S), 1) ;
%     plot(t, YS.lower, 'r') ;
%     plot(t, YS.upper, 'r') ;
% end
% 
% if ~isempty(S) || ~isempty(E)
%     YE = join(drugs.compute([S ; E]), 1) ;
%     plot(t, YE.lower, 'g') ;
%     plot(t, YE.upper, 'g') ;
% end
% 
% draw_boxes(S,E,N,1,2) 


%% Solving an IVGTT problem (High complexity: more than 1h computations time and 10 Gb RAM memory required)
% In this fourth and last example, an intravenous glucose tolerance test
% (IVGTT) is considered. Thorough explanations about it are given in
% <../examples/html/IVGTT_parameters.html IVGTT_parameters> but, to put it
% simply, it is the same kind of problem than the previous application,
% that is a 4-D data-fitting problem, except that it is slightly trickier,
% because:
% 
% * No explicit expression of the function to be inverted is known;
% integrating this function is required in order to determine its value in
% a given point
%
% * The theoretical model is not accurate enough to fit all the
% measurements, especially at the beginning of the test
%
% As a consequence, in <../examples/html/IVGTT_parameters.html IVGTT_parameters>,
% a numerical integration is performed in the _compute_ function, that
% involves therefore much more computations than in the previous
% applications. In addition, modal interval analysis (MIA) has been used in
% order to reduce overbounding phenomena in interval computations. At last,
% the (optional) parameter _tol___out_ in <vsivia_parameters.html
% vsivia_parameters>, that stands for the number of measurements that may
% not be matched, has been set to 6.
%
% As in the previous application about a bolus intravenous injection of a
% drug into a subject, a graphical representation of the model has been
% drawn using the enclosures of the parameters computed by <vsivia.html
% vsivia>.

% ivgtt = IVGTT_parameters ;
% 
% [S,E] = vsivia(ivgtt) ;
% 
% figure ;
% 
% hold on ;
% 
% t = ivgtt.t ;
% 
% Y0 = ivgtt.Y0 ;
% 
% errorbar(t, .5*(Y0.lower+Y0.upper), .5*(Y0.upper-Y0.lower), 'kx') ;
% 
% YS = join(ivgtt.compute(S), 1) ;
% 
% plot(t, YS.lower, 'r') ;
% plot(t, YS.upper, 'r') ;
% 
% YE = join(ivgtt.compute([S ; E]), 1) ;
% 
% plot(t, YE.lower, 'g') ;
% plot(t, YE.upper, 'g') ;

