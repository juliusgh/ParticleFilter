function [S, E, N, Y0, C] = vsivia(params)

S = interval([]) ;
E = interval([]) ;
N = interval([]) ;
C = interval([]) ;

if ~isa(params, 'vsivia_parameters')
   error('Wrong argument ''params'' (must be an instance of vsivia_parameters)' ) ;
end


if strcmp(params.algorithm, 'inversion')
    inversion = true ;
    optimization = false ;
    use_local_opt = false ;
elseif strcmp(params.algorithm, 'optimization')
    inversion = false ;
    optimization = true ;
    use_local_opt = false ;
elseif strcmp(params.algorithm, 'fast optimization')
    inversion = false ;
    optimization = true ;
    use_local_opt = true ;
else
    error('Wrong argument ''algorithm'' (must be either ''inversion'' or ''optimization'' or ''fast optimization'')') ;
end


U0 = params.U0 ;

if isa(U0, 'interval')
    U = U0 ;
elseif isa(U0, 'numeric')
    U = interval(U0, [], 2)' ;
else
    error('Wrong argument ''U0'' (must be an interval or an array of reals)') ;
end

w0 = width(U(1,:)) ;

[epsilon, absolute, relative, volume] = getEpsilonParams(params.epsilon) ;

maxBoxes = params.maxBoxes;

epsbisect = absolute | relative ;
epsresult = volume ;

epsilon(isnan(epsilon)) = inf ;

if ~all(epsilon >= 0)
    error('Wrong argument ''epsilon'' (must be positive)') ;
end


contract_kind = 'none' ;

try
    contract_kind = params.contraction_kind ;
catch err
end


contract_thres = .5 ;

try
    contract_thres = params.contraction_threshold ;
catch err
end


if ~isa(contract_thres, 'numeric') || ~(contract_thres >= 0)
    error('Wrong argument ''contraction_threshold'' (must be a positive real)') ;
end

if isa(contract_kind, 'char')
    
    if strcmp(contract_kind, 'none')
        
        use_contraction = false ;
        
    elseif strcmp(contract_kind, 'auto')
        
        use_contraction = true ;
        
        contractFun = str2func(writeFW(params)) ;

%     elseif strcmp(contract_kind, 'custom')
%         
%         if sum(strcmp(methods(params), 'contract')) > 0
%         
%             use_contraction = true ;
% 
%             contractFun = @params.contract ;
%             
%             nincontract = nargin(str2func([class(params) '.constraint'])) ;
%     
%             contract_box = (nincontract == 1) | (nincontract == -2) ;
%             
%         else
%             
%             error('Custom contraction specified but ''contract'' function missing...') ;
%         
%         end
        
    else
        
        error('Wrong argument ''contraction_kind'' (must be either ''none'' or ''auto'')') ;
        
    end
    
else
    
    error('Wrong argument ''contraction_kind'' (must be either ''none'' or ''auto'' or ''custom'')') ;
    
end


Y0 = interval.toInterval(params.Y0) ;
     
if use_contraction
    
     nIn = abs(nargin(contractFun)) ;
     nOut = abs(nargout(contractFun)) ;
     
     if abs(nIn) < size(U,2) + 1
         
         disp('Wrong contraction function (incorrect number of arguments)') ;
         
         use_contraction = false ;
         
     end
     
end

nincompute = length(params.getParams('compute')) ;

compute_box = (nincompute == 1) | (nincompute == -2) ;

if use_local_opt && ~compute_box
    error('Fast optimization requires a ''compute'' method whose input is a whole box') ;
end


if sum(strcmp(methods(params), 'constraint')) > 0
    nonlin_constraints = @params.constraint ;
else
    nonlin_constraints = [] ;
end


if ~isempty(nonlin_constraints)
    
    ninconstraint = nargin(str2func([class(params) '.constraint'])) ;
    
    constraint_box = (ninconstraint == 1) | (ninconstraint == -2) ;
    
    if use_local_opt && ~constraint_box
        
        error('Fast optimization requires a constraint function whose input is a whole box') ;
        
    end
    
end


A = [] ;

try
    A   = params.A ;
catch err
end

A = [] ;

try
    A   = params.A ;
catch err
end

Aeq = [] ;

try
    Aeq   = params.Aeq ;
catch err
end

b = [] ;

try
    b   = params.b ;
catch err
end

beq = [] ;

try
    beq   = params.beq ;
catch err
end


Amat = interval_matrix2D(A, [], 0) ;
Aeqmat = interval_matrix2D(Aeq, [], 0) ;


constrained = true ;

try
    constrained = params.constrained ;
catch err
end

if ~isa(constrained, 'logical')
   error('Wrong argument ''constrained'' (must be either true or false)') ; 
end


lb = interval(U.lower, [], 0) ;
ub = interval(U.upper, [], 0) ;


if size(A,1) ~= size(b,1)
    error('Dimensions of A and b inconsistent with each other') ;
end

if size(Aeq,1) ~= size(beq,1)
    error('Dimensions of A and b inconsistent with each other') ;
end

if constrained && isempty(A) && isempty(b) && isempty(Aeq) && isempty(beq)...
        && isempty(lb) && isempty(ub) && isempty(nonlin_constraints)
    
    constrained = false ;
    
end


tol_min = 1e-8 ;

try
    if isa(params.tol_min, 'numeric') && params.tol_min >= 0
        tol_min = params.tol_min ;
    else
        error('Wrong argument ''tol_min'' (must be a positive real number)') ;
    end
catch err
end


tol_in = 0 ;

tol_in_abs = true ;
tol_in_rel = false ;

try
    [tol_in, tol_in_abs, tol_in_rel] = getTolParams(params.tol_in, 'tol_in') ;
catch err
end


tol_out = 0 ;

tol_out_abs = true ;
tol_out_rel = false ;

try
    [tol_out, tol_out_abs, tol_out_rel] = getTolParams(params.tol_out, 'tol_out') ;
catch err
end


VS = 0 ;

nbox = 0 ;

nit = 0 ;

perform_fmin = use_local_opt ;

oldmin = inf ;

np = size(U,2) ;


st = tic ;

while size(U,1) > 0
    
    nit = nit + 1 ;
    if nit > 25
        disp(['iteration: ',num2str(nit)]);
    end
    
    nbox = nbox + size(U,1) ;
    
    if constrained
        
        out = false(size(U,1), 1) ;
        
        if ~isempty(Amat)
            out(Amat*U' > b) = true ;
        end
        
        if ~isempty(Aeqmat)
            out(isout(beq,Aeqmat*U')) = true ;
        end
        
        if ~isempty(lb)
           out(sum(U < lb, 2) > 0) = true ;
        end
        
        if ~isempty(ub)
           out(sum(U > ub, 2) > 0) = true ;
        end
        
        if ~isempty(nonlin_constraints)
            
            if constraint_box
                [c, ceq] = params.constraint(U) ;
            else
                compute_in = cell(1, size(U,2)) ;
                [compute_in{:}] = cell(U,1) ;
                [c, ceq] = params.constraint(compute_in{:}) ;
            end
            
            out(any(c > 0, 2)) = true ;
            out(any(ceq ~= 0, 2)) = true ;
            
        end
        
        C = cat(1, C, U(out,:)) ;
        
        U = U(~out,:) ;
        
    end
    
    if compute_box
        Y = params.compute(U) ;
    else
        compute_in = cell(1, size(U,2)) ;
        [compute_in{:}] = cell(U,1) ;
        Y = params.compute(compute_in{:}) ;
    end
    
    if inversion
        
        if tol_in_rel
            in = sum(~isin(Y, Y0),2)/size(Y,2) <= tol_in ;
        elseif tol_in_abs
            in = sum(~isin(Y, Y0),2) <= tol_in ;
        end
            
        if tol_out_rel
            out = sum(isout(Y, Y0),2)/size(Y,2) > tol_out ;
        elseif tol_out_abs
            out = sum(isout(Y, Y0),2) > tol_out ;
        end
        
    elseif optimization
        
        if isa(Y, 'interval')
            minlow = min(Y.lower, [], 1) ;
            minupp = min(Y.upper, [], 1) ;
        else
            minlow = min(Y, [], 1) ;
            minupp = minlow ;
        end
        
        Y0 = Y0 & interval(minlow - eps(minlow), minupp + eps(minupp)) ;
        
        in = false(size(U,1), 1) ;
        
        out = isout(Y,Y0) ;
        
    end
    
    if epsbisect
        
        if size(S,1) + size(E,1) + size(U,1) > maxBoxes

            epsilon = Inf;

        end
        
        S = cat(1, S, U(in,:)) ;
        
        if nargout > 2
            N = cat(1, N, U(out,:)) ;
        end
        
        U = U(~(in|out),:) ;
        
        tiny = min(bsxfun(@lt, width(U), epsilon), [], 2) ;
        
        E = cat(1, E, U(tiny,:)) ;
        
        U = U(~tiny,:) ;
        
    end
    
    if epsresult
        
        DS = U(in,:) ;
        
        VS = VS + sum(volume(DS,2)) ;
        
        S = cat(1, S, DS) ;
        
        if nargout > 2
           N = cat(1, N, U(out,:)) ; 
        end
        
        U = U(~(in|out),:) ;
        
        VU = sum(volume(U,2)) ;
        
        if VU/VS < epsilon
            
            E = U ;
            
            break ;
                        
        end
        
    end
        
    n = size(U,1) ;
    
    if n > 0
        
        if perform_fmin
            
            wmin = min(epsilon(:)) ;
            
            options = optimset('Display', 'off', 'LargeScale', 'off', 'Algorithm', 'sqp', 'TolX', wmin) ;
            
            if constrained
                
                xmin = cell2mat(cellfun(@getMinCon, num2cell(mid(U),2), 'UniformOutput', false)) ;
                
            else
                
                xmin = cell2mat(cellfun(@getMinUnc, num2cell(mid(U),2), 'UniformOutput', false)) ;
                
            end
            
            xmin = sortrows(xmin, 1) ;
            
            minval = min(oldmin,xmin(1,1)) ;
            
            if abs((oldmin-minval)/oldmin) < tol_min
                
                perform_fmin = false ;
                
            else
                
                oldmin  = minval ;
                
                xmin = xmin(1, 2:end) ;
                
                Y0 = Y0 & interval([-inf params.compute(xmin)]) ;
                
            end
            
        end
        
        
        if use_contraction
            
            inArgs = cell(1,nOut) ;
            
            if compute_box
                inArgs{1} = U ;
                inArgs(2:nOut) = {interval([-inf inf])} ;
            else
                [inArgs{1:np}] = cell(U,1) ;
                inArgs(np+1:nOut) = {interval([-inf inf])} ;
            end
            
            inArgs{end} = Y0 ;
            
            outArgs = cell(1,nOut) ;
            
            s = true(n,1) ;
            
            while sum(s) > 0
                
                [outArgs{:}] = feval(contractFun, inArgs{:}) ;
                
                if compute_box
                    boxout = outArgs{1} ;
                    wout = width(boxout) ;
                    xcontract = (min(wout ./ width(inArgs{1}), [], 2) < contract_thres)...
                                    & ~min(bsxfun(@lt, wout, epsilon), [], 2) ;
                else
                    boxout = horzcat(outArgs{1:np}) ;
                    wout = width(boxout) ;
                    xcontract = (min(wout ./ width(horzcat(inArgs{1:np})), [], 2) < contract_thres)...
                                    & ~min(bsxfun(@lt, wout, epsilon), [], 2) ;
                end
                                                                         % xcontract: Continue contraction, among selected
                nxcontract = ~xcontract ;                                % nxcontract: Do not continue contraction, among selected
                
                contract = zeros(n,1) ;
                
                contract(s) = xcontract ;                                % contract : continue contraction, among all
                
                ncontract = ~contract ;                                  % ncontract : do not continue contraction, among all
                
                push = s & ncontract ;                                   % push: do not continue contraction for selected only, among all
                
                if sum(push) > 0
                    U(push,:) = boxout(nxcontract,:) ;
                end
                
                inArgs = cellfun(@(x) x(xcontract,:) , outArgs, 'UniformOutput', false) ;
                
                s = s & contract ;
                
            end
            
            U = U(sum(U.lower > U.upper,2) == 0, :) ;
            
        end
        
    end
    
    if ~isempty(U)
        
        U = bisect(U,w0,2,1) ;
        
    end
            
end


if optimization && constrained
    
    in = true(size(E,1), 1) ;
    
    if ~isempty(Amat)
        in(~(Amat*E' <= b)) = false ;
    end
    
    if ~isempty(Aeqmat)
        in(~isin(beq,Aeqmat*E')) = false ;
    end
    
    if ~isempty(lb)
        in(~all(E >= lb, 2)) = false ;
    end
    
    if ~isempty(ub)
        in(~all(E <= ub, 2)) = false ;
    end
    
    if ~isempty(nonlin_constraints)
        
        if constraint_box
            [c, ceq] = params.constraint(E) ;
        else
            compute_in = cell(1, size(E,2)) ;
            [compute_in{:}] = cell(E,1) ;
            [c, ceq] = params.constraint(compute_in{:}) ;
        end
        
        in(~all(c <= 0, 2)) = false ;
        in(~all(ceq == 0, 2)) = false ;
        
    end
    
    S = E(in,:) ;
    
    E = E(~in,:) ;
    
end
    

toc(st) ;

disp(['Number of iterations: ' num2str(nit)]) ;
    
disp(['Number of boxes processed: ' num2str(nbox)]) ;


    function y = getMinUnc(x)
        
        [x, val] = fminunc(@params.compute, x, options) ;
        
        y = [val, x] ;
        
    end


    function y = getMinCon(x)
        
        [x, val] = fmincon(@params.compute, x, A, b, Aeq, beq, lb.lower, ub.lower, nonlin_constraints, options) ;
        
        y = [val, x] ;
        
    end


    function [epsilon, absolute, relative, volume] = getEpsilonParams(epsilon_param)
        
        if isa(epsilon_param, 'numeric')
            
            epsilon = epsilon_param ;
            
            absolute = true ;
            relative = false ;
            volume   = false ;
            
        elseif isa(epsilon_param, 'char')
            
            parts = regexp(epsilon_param, '(?<number>.*?)(?<suffix>%|rel|vol)?\s*$', 'tokens') ;
            
            parts = parts{:} ;
            
            epsilon = str2num(parts{1}) ;         %#ok<ST2NM> - parsing of arrays required
            
            if numel(epsilon) == 0
                
                error('Wrong argument ''epsilon'' (must contain a number or an array of numbers parsable by ''str2num'')') ;
                
            end
            
            kind   = parts{2} ;
            
            absolute = strcmp(kind, '') ;
            relative = strcmp(kind, '%') | strcmp(kind, 'rel') ;
            volume   = strcmp(kind, 'vol') ;
            
            if strcmp(kind, '%')
                epsilon = .01 * epsilon ;
            end
            
        else
            
            error('Wrong argument ''epsilon'' (type must be number or char)') ;
            
        end
        
    end


    function [tol, absolute, relative] = getTolParams(tol_param, tol_kind)
        
        if isa(tol_param, 'numeric')
            
            tol = tol_param ;
            
            absolute = true ;
            relative = false ;
            
        elseif isa(tol_param, 'char')
            
            parts = regexp(epsilon_param, '(?<number>.*?)(?<suffix>%|rel)?\s*$', 'tokens') ;
            
            parts = parts{:} ;
            
            tol = str2double(parts{1}) ;
            
            if isnan(tol)
                
                error(['Wrong argument ''' tol_kind ''' (must contain a valid number parsable by ''str2double'')']) ;
                
            end
            
            kind   = parts{2} ;
            
            absolute = strcmp(kind, '') ;
            relative = strcmp(kind, '%') | strcmp(kind, 'rel') ;
            
            if strcmp(kind, '%')
                tol = .01 * tol ;
            end
            
        else
            
            error(['Wrong argument ''' tol_kind ''' (type must be number or char)']) ;
            
        end
        
    end


end
