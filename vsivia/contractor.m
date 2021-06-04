%% Generating automatically a forward-backward contractor
% The class _contractor_, presented hereafter, is used by the main
% algorithm, _vsivia_, in order to generate in an automatic way a
% forward-backward contractor from a Matlab function using only features
% (_i.e._ methods and operators) defined by this class (_contractor_).

%% Approach
% The process relies upon operator and function overloading.
% Originally, the Matlab function to be processed is intented to work
% with scalars and/or intervals (see class _interval_). When executed,
% the Matlab interpreter calls the methods specified in this function,
% one after another. In the case of scalars or intervals, it results
% in a sequence of artihmetic operations. Here, the idea consists in
% using those method calls in order to build the forward-backward contractor,
% step-by-step, instead of performing these arithmetic operations.
%
% Concretly, it means that, at each step (_i.e._ each operation), the
% intructions corresponding to the forward and the backward are generated,
% and then put together with the intructions previously generated. Those
% list of instructions are built separately. The forward list is a FIFO list,
% while the backward one is a LIFO one. Moving back on the building process,
% each one of its steps is given a number, that is used to define intermediate
% variables in a unique way. In parallel, an other list containing all the
% intermediate variables beginning by a "A" (as opposed to those beginning by a
% "B") introduced so far is built. This list may be used to generate a header
% for the function designed to contain the forward-backward contractor. 
%
% About the names of the intermediate variables, those beginning with a "A"
% correspond to the result of a step of the forward process, while those starting
% with a "B" are used internally, in order to store the value of a particular variable
% (an input variable or a type "A" intermediate variable) during the forward
% process, before it is applied a non-injective operation (e.g. the square function
% or a trigonometric function). During the backward process, "B" variables help,
% in some cases, to perform the corresponding inverse operations with a more accurate
% result.


classdef contractor
    
    properties
        
        name,                    % The step number
        
        forward = {},            % List of the forward instructions
        
        backward = {},           % List of the backward instructions
        
        intermediate = {} ;      % Intermediate variables list
        
    end
    
    
    %% Giving intermediate variables unique names
    % Any intermediate variable is given a name of the form "A" or "B",
    % followed by a stricly positive integer (_e.g._ "A1" or "B12").
    % Numbers related to type "A" or "B" intermediate variables are generated
    % using the two methods below, that basically increment an internal counter
    % each time they are called, and return a name built using the value of this counter.
    
    methods (Hidden, Static)
        
        function s = nextA()                   % Generates a name for type "A" intermediate variables
            
            persistent n ;                     % n : internal counter
            
            if isempty(n)                      % Initializes n to 1 if empty (first call)
                n = 1 ;
                s = 'A1' ;
            else
                n = n+1 ;                      % Otherwise, increments it and builds the name
                s = ['A' num2str(n, '%i')] ;
            end
            
        end
        
        function s = nextB()                   % Generates a name for type "A" intermediate variables
            
            persistent n ;                     % n : internal counter
            
            if isempty(n)                      % Initializes n to 1 if empty (first call)
                n = 1 ;
                s = 'B1' ;
            else
                n = n+1 ;                      % Otherwise, increments it and builds the name
                s = ['B' num2str(n, '%i')] ;
            end
            
        end

        
    end
    
    
    %%
    % -- *Miscellaneous* --
    %
    % Some binary methods working with _contractor_ instances (_e.g._ _plus_, _minus_)
    % also accept non-_contractor_ parameters, for instance numbers. At the beginning
    % of such methods, input parameters are ensured to be _contractor_ instances
    % (and if not, transformed to) by calling the following method.

    methods (Static)
        
        function c = toContractor(c)
           
            if ~isa(c, 'contractor')
                c = contractor(c) ;
            end
            
        end
        
    end % static methods
    
    
    %% Instantiating _contractor_ objects
    %
    % The constructor hereafter can be invoked using different combinaisons
    % of parameters. In any case, a first parameter is required and shall be:
    % * A _contractor_ instance: acts as copy constructor
    % * A number: creates a _contractor_ object whose name will be the string representation of the number
    % * A string: creates a _constractor_ object whose name will be that string
    % Then, a second and even a third parameter can be specified. They are optional, shall be
    % instances of _contractor_, and are used to built the _intermediate_ property of the
    % newly created _contractor_. The _intermediate_ property of the "child" contains
    % every unique value contained in that (those) of its parent(s).
    
    
    methods
        
        function c = contractor(varargin)
            
            error(nargchk(1,3,nargin)) ;
            
            arg = varargin{1} ;
            
            if isa(arg, 'contractor')
                
                c = arg ;
                
            elseif isa(arg, 'numeric')
                
                if length(arg) > 1
                
                    c.name = ['[ ' num2str(arg, '%20.20g') ' ]'] ;
                    
                else
                    
                    c.name = num2str(arg, '%20.20g') ;
                    
                end
                
            else
                
                c.name = arg ;
                
            end
            
            if nargin == 2
                
                if c.name(1) == 'A'
                
                    c.intermediate = [ varargin{2}.intermediate ; { c.name } ] ;
                    
                end
                
            elseif nargin == 3
                
                c.intermediate = [ unique([ varargin{2}.intermediate ;
                                              varargin{3}.intermediate ]) ;
                                   { c.name } ] ;
                               
            end
            
        end % contractor
        
        
        %% Generating forward and backward instructions in a generic way
        % The methods below are intented to generate the forward and backward
        % properties of a "child" from that (those) of its "parent(s)" (_i.e._,
        % when a operation is performed, the "parent(s)" is (are) the operand(s),
        % while the "child" is the result. Four cases are considered, depending
        % upon the number of operands (one or two) and the forward or backward
        % nature of the property to be built.
        
        function s = genBinaryForward(r, a, b, op)
            s = [ a.forward ;
                  b.forward ;
                  { [ r.name ' = ' r.name ' & (' a.name op b.name ') ;' ] } ] ;
        end % genBinaryForward
        
        function s = genBinaryBackward(r, a, b, op1, op2, inv)
            
            s = {} ;
            
            if isvarname(a.name)
                s = { [ a.name ' = ' a.name ' & (' r.name op1 b.name ') ;' ] } ;
            end
            
            if isvarname(b.name)
                if inv
                    s = [ s ; { [ b.name ' = ' b.name ' & (' a.name op2 r.name ') ;' ] } ] ;
                else
                    s = [ s ; { [ b.name ' = ' b.name ' & (' r.name op2 a.name ') ;' ] } ] ;
                end
            end
            
            s = [ s ; a.backward ; b.backward ; ] ;
            
        end % genBinaryBackward
        
        function s = genUnaryForward(r, a, fun)
            s = [ a.forward ;
                  { [ r.name ' = ' r.name ' & ' fun '(' a.name ') ;' ] } ] ;
        end % genUnaryForward
        
        function s = genUnaryBackward(r, a, fun)
            s = [ { [ a.name ' = ' a.name ' & ' fun '(' r.name ') ;' ] } ;
                  a.backward ; ] ;
        end % genUnaryBackward
        
        
        function r = plus(a, b)
            a = contractor.toContractor(a) ;
            b = contractor.toContractor(b) ;
            r = contractor(contractor.nextA, a, b) ;
            r.forward = genBinaryForward(r, a, b, '+') ;
            r.backward = genBinaryBackward(r, a, b, '-', '-',false) ;
        end % plus
        
        function r = minus(a,b)
            a = contractor.toContractor(a) ;
            b = contractor.toContractor(b) ;
            r = contractor(contractor.nextA, a, b) ;
            r.forward = genBinaryForward(r, a, b, '-') ;
            r.backward = genBinaryBackward(r, a, b, '+', '-',true) ;
        end % minus
        
        function r = uminus(a)
            r = contractor(contractor.nextA, a) ;
            r.forward =  [ a.forward ;
                           { [ r.name ' = ' r.name ' & -' a.name ' ;' ] } ] ;
            r.backward = [ { [ a.name ' = ' a.name ' & -' r.name ' ;' ] } ;
                           a.backward ] ;
        end % uminus
        
        function r = mtimes(a,b)
            a = contractor.toContractor(a) ;
            b = contractor.toContractor(b) ;
            r = contractor(contractor.nextA, a, b) ;
            r.forward = genBinaryForward(r, a, b, '*') ;
            r.backward = genBinaryBackward(r, a, b, '/', '/',false) ;
        end % mtimes
        
        function r = mrdivide(a,b)
            a = contractor.toContractor(a) ;
            b = contractor.toContractor(b) ;
            r = contractor(contractor.nextA, a, b) ;
            r.forward = genBinaryForward(r, a, b, '/') ;
            r.backward = genBinaryBackward(r, a, b, '*', '/',true) ;
        end % mrdivide
        
        function r = exp(a)
            r = contractor(contractor.nextA, a) ;
            r.forward = genUnaryForward(r, a, 'exp') ;
            r.backward = genUnaryBackward(r, a, 'log') ;
        end % exp
        
        function r = cos(a)
            
            r = contractor(contractor.nextA, a) ;
            
            r.forward = genUnaryForward(r, a, 'cos') ;
            
            bname = contractor.nextB ;
            
            r.forward = [ r.forward ; { [ bname ' = ' a.name ' ;' ] } ] ;
            
            r.backward = [ { [ a.name ' = ' a.name ' & acos(' r.name ',' bname ') ;' ] } ;
                             a.backward ] ;
            
        end % cos
        
        function r = sin(a)
            
            r = contractor(contractor.nextA, a) ;
            
            r.forward = genUnaryForward(r, a, 'sin') ;
            
            bname = contractor.nextB ;
            
            r.forward = [ r.forward ; { [ bname ' = ' a.name ' ;' ] } ] ;
            
            r.backward = [ { [ a.name ' = ' a.name ' & asin(' r.name ',' bname ') ;' ] } ;
                             a.backward ] ;
                         
        end % sin
        
        
        function r = mpower(a,n)
            
            r = contractor(contractor.nextA, a) ;
            
            r.forward = [ a.forward ;
                          { [ r.name ' = ' r.name ' & ' a.name '^' num2str(n, '%20.20g') ' ;' ] } ] ;
                      
            if n == round(n)
                
                if mod(n,2) == 0
                    
                    bname = contractor.nextB ;
                    
                    r.forward = [ r.forward ; { [ bname ' = ' a.name ' ;' ] } ] ;
                    
                    r.backward = [ { [ a.name ' = ' a.name ' & nthroot(' r.name ','...
                                        num2str(n, '%20.20g') ',' bname ') ;' ] } ;
                                        a.backward ; ] ;
                                 
                else
                    
                    r.backward = [ { [ a.name ' = ' a.name ' & nthroot(' r.name ',' num2str(n, '%20.20g') ') ;' ] } ;
                                        a.backward ; ] ;
                    
                end
                
            else
                
                r.backward = [ { [ a.name ' = ' a.name ' & ' r.name '^' num2str(1/n, '%20.20g') ' ;' ] } ;
                                 a.backward ; ] ;
                
            end
            
        end % mpower
        
        function r = power(r,n)
            r = mpower(r,n) ;
        end % power
        
        
        function r = subsref(a,s)
            
            s = s(1) ;
            
            if strcmp(s.type, '.')
                
                 r = builtin('subsref', a, s) ;
                 
            elseif strcmp(s.type, '()')
                
                r = contractor(contractor.nextB, a) ;
                
                r.forward = [ a.forward ;
                    { [ r.name ' = subsref(' a.name ', cell2struct({'...
                    '''' s.type ''' ; ' subs2str(s.subs)...
                    '}, {''type'', ''subs''}, 1)) ;' ] } ] ;
                
                r.backward = [ { [ a.name ' = subsasgn(' a.name ', cell2struct({'...
                    '''' s.type ''' ; ' subs2str(s.subs)...
                    '}, {''type'', ''subs''}, 1), ' r.name ') ;' ] } ;
                    a.backward ] ;
                
            end
            
        end % subsref
        
        
    end % methods
    
end

