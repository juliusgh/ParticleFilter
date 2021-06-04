%% Specifying parameters for _vsivia_

%% Overview
% The function <vsivia.html vsivia> takes, as its input, an object
% implementing this abstract class. Its methods and properties give
% <vsivia.html vsivia> a complete description of the problem to be solved. 


classdef vsivia_parameters
    
    %% Mandatory parameters
    % The following parameters define the problem itself; they shall
    % therefore be defined by the user. They are declared as "abstract" for
    % this purpose.
    
    properties (Abstract)
        
        %%
        % *algorithm* : tells <vsivia.html vsivia> about the kind of
        % problem and the algorithm to be used in order to solve it:
        %
        % * _inversion_ : for any inversion problem
        %
        % * _optimization_ : for any optimization problem
        %
        % * _fast optimization_ : also for optimization problems; uses
        % local optimization in order to find an upper bound to the global
        % minimum more quickly (this may not work efficiently in any case)
        
        algorithm ;
        
        %%
        % *U0* : specifies the initial box, that is a 1-by-n interval or an
        % n-by-2 array of doubles
        
        U0 ;
        
        %%
        % *Y0* : specifies the set to be inverted (for an inversion
        % problem) or an initial guess of the global minimum (_i.e._
        % an interval containing it, for optimization problems).
        
        Y0 ;
        
        %%
        % *epsilon* : the accuracy parameter to be used, that tells
        % <vsivia.html vsivia> when to stop bisecting boxes;
        % may be either:
        %
        % * A number (integer or double): defines a minimal size below
        % which a box cannot be bisected any more
        %
        % * A number (integer or double) within a string: same thing as
        % for the previous point
        %
        % * A number followed by "rel" or "%", within a string: defines the
        % minimal size with respect to the size of _U0_ ; "rel" indicates a
        % scale factor, while "%" indicates a percentage (thus "0.04 rel"
        % is equivalent to "4 %").
        %
        % * A number followed by "vol", within a string: considers the
        % ratio between the volume of the current frontier (_i.e._ boxes
        % that shall be bisected) and that of the current solution set
        % (_i.e._ boxes that have been proven to be solutions), rather than
        % the size of the boxes that should be bisected; iterations are
        % stopped when that ratio becomes smaller than the specified
        % number
        %
        % Please note that:
        %
        % * consecutive spaces count like a single space
        %
        % * _epsilon_ can be either a single value, or a vector of the size
        % of _U0_ ; in this last case, each component of the problem has a
        % dedicated threshold
        %
        % Multiple values, for a vector, shall be specified in a same
        % string. The kind of accuracy parameter (absolute, "rel", "%" or
        % "vol") is the same for every component of the problem, whether
        % _epsilon_ be a single value or a vector. There is at most one of
        % those suffixes, placed at the end of the string.
        
        epsilon ;
        
    end % abstract properties
    
    
    %% Optional parameters
    % The following parameters can be defined, as properties, in order to
    % tune the solving algorithm in a finer way:
    %
    % * contraction_kind
    %
    % * contraction_threshold
    %
    % * tol_min (full vsivia)
    %
    % * tol_in
    %
    % * tol_out
    %
    % * constrained (full vsivia)
    %
    % * A (full vsivia)
    %
    % * Aeq (full vsivia)
    %
    % * b (full vsivia)
    %
    % * beq (full vsivia)
    
    methods (Abstract)
        
        y = compute(varargin) ;
        
    end
    
    
    methods
        
        function output = genFW(vsivia_params)
            
            r = getParams(vsivia_params, 'compute') ;
            
            params = cellfun(@(s) contractor(s), r, 'UniformOutput', false) ;
            
            r = vsivia_params.compute(params{:}) ;
            
            params = cellfun(@(x) x.name, params', 'UniformOutput', false) ;
            
            vars = [ params ; r.intermediate ] ;
            
            vars(1:end-1) = cellfun(@(s) [s ', '], vars(1:end-1), 'UniformOutput', false) ;
            
            vars = cell2mat(vars') ;
            
            output{1,1} = [ 'function [' vars '] = auto_forward_backward(' vars ', varargin)' ] ;
            output{2,1} = '' ;
            output{3,1} = '% Forward' ;
            output{4,1} = '' ;
            
            output = [ output ; r.forward ] ;
            
            if strcmp(vsivia_params.algorithm, 'optimization')
                output{end+1,1} = '' ;
                output{end+1,1} = ['empty = ~(' r.name '.lower <= min(' r.name '.upper, [], 1)) ;' ] ;
                output{end+1,1} = [ r.name '.lower(empty) = inf ;' ] ;
                output{end+1,1} = [ r.name '.upper(empty) = -inf ;' ] ;
                output{end+1,1} = '' ;
            end
            
            output{end+1,1} = '' ;
            output{end+1,1} = '% Backward' ;
            output{end+1,1} = '' ;
            
            output = [ output ; r.backward ] ;
            
            output{end+1,1} = '' ;
            output{end+1,1} = 'end' ;
            
        end
        
        
        function [funname, fullpath] = writeFW(vsivia_params)
            
            funname = 'auto_forward_backward' ;
            
            fullpath = strcat(fileparts(mfilename('fullpath')), '/auto_forward_backward.m') ;
            
            f = fopen(fullpath, 'w+') ;
            
            output = genFW(vsivia_params) ;
            
            fprintf(f, '%s\r\n', output{:}) ;
            
            fclose(f) ;
            
            % Refreshes newly created file for Matlab
            
            clear auto_forward_backward ;
            
            exist('auto_forward_backward.m', 'file') ;
            
        end
        
        
        function r = getParams(vsivia_params, funName)
            
            m = methods(vsivia_params, '-full') ;
            
            r = cellfun(@(s) regexp(s, '\s|\(|\)|,', 'split'), m, 'UniformOutput', false) ;
            
            static = cellfun(@(s) strcmp(s{1}, 'Static'), r) ;
            
            r(static) = cellfun(@(s) s(2:end), r(static), 'UniformOutput', false) ;
            
            fun = cellfun(@(s) s{2}, r, 'UniformOutput', false) ;
            
            r = r{strcmp(fun, funName)} ;
            
            r = r(3:end) ;
            
            r = r(cellfun(@(s) ~isempty(s), r)) ;
            
            if ~static(strcmp(fun, funName))
                
                r = r(2:end) ;
                
            end
            
        end
        
    end
    
end
