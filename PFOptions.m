classdef PFOptions
    %PFOPTIONS Options for Particle Filters
    %   Detailed explanation goes here
    
    properties
        samples % particle count
        iterations % iterations of filtering
        simulate % simulate xvals, zvals (logical)
        xvals % fixed system state values
        zvals % fixed measurement values
        eps % edge-length of grid
        debug % debug output (logical)
        LPOptions % options for linprog
        gurobi % use gurobi linprog solver (logical)
        reduceToHull % store only convex hull instead of all particles (logical)
    end
    
    methods
        function obj = PFOptions(samples,iterations,simulate,xvals,zvals)
            %PFOPTIONS Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 3
                simulate = true;
                xvals = [];
                zvals = [];
                if nargin < 2
                    iterations = 16;
                    if nargin < 1
                        samples = 5000;
                    end
                end
            end
            obj.samples = samples;
            obj.iterations = iterations;
            obj.simulate = simulate;
            obj.xvals = xvals;
            obj.zvals = zvals;
            obj.debug = true;
            obj.LPOptions = optimoptions('linprog','Algorithm','dual-simplex',...
                    'OptimalityTolerance',1e-7,'ConstraintTolerance',1e-4,...
                    'Preprocess','basic','Display','none');
            obj.gurobi = true;
            obj.reduceToHull = false;
        end
    end
end

