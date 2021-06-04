%PS_VSIVIA Set-Membership Filtering via Interval Analysis
%   for 4-dim. ManipulatorSystem/PendulumSystem

maxk = 51; % iterations
init = true;
if init
    S = cell(1,maxk); % solution boxes
    E = cell(1,maxk); % undefined boxes
    N = cell(1,maxk); % non-solution boxes
    pfv = cell(1,maxk); % vsivia parameters
    Y = cell(1,maxk); % measurement intervals
end
P = []; % propagated boxes from last iteration
intersect = false; % intersection mode (logical)
U0 = system.X0.bounds; % start interval for vsivia
V = system.V.bounds; % process noise interval
u = mpc1R.u; % hier Stellgrößen der Regelung laden!
%u = zeros(2,1);
propagate = @(x,k) (system.gpredict(x,u(:,k))' + V)'; % function handle for propagation
measure = @system.h; % function handle for measurement

for k = startk:maxk
   fprintf('k = %d\n',k);
   Y{k} = system.z(:,k) + system.W.bounds;
   pfv{k} = PF_vsivia(measure,Y{k},U0,P,intersect);
   pfv{k}.epsilon = '2%';
   if k == startk
       pfv{k}.epsilon = '1%';
   end
   [S{k},E{k},N{k}] = vsivia(pfv{k});
   X = [S{k} E{k}];
%    if k == startk
%        X = [X E{k}];
%    end
   P = propagate(X,k); % propagate solution boxes
   U0 = P.box() + 0.05*interval(-ones(1,2),ones(1,2)); % set new start interval for vsivia
   intersect = true; % activate intersection mode for next iteration
end