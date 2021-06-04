%PS_VSIVIA Set-Membership Filtering via Interval Analysis
%   (only generates set H_k)

maxk = 10; % iterations
Sz = cell(1,maxk); % solution boxes
Ez = cell(1,maxk); % undefined boxes
Nz = cell(1,maxk); % non-solution boxes
psv = cell(1,maxk); % vsivia parameters
Y = cell(1,maxk); % measurement intervals
% P = []; % propagated boxes from last iteration
% intersect = false; % intersection mode (logical)
U0 = system.x0(1:2) + interval([-0.1 -0.1],[0.1 0.1]); % start interval for vsivia
% padding = interval([-0.1 -0.1],[0.1 0.1]); % padding around start interval
% V = system.V.toInterval(); % process noise interval
% propagate = @(x) system.g(x) + V; % function handle for propagation
measure = @system.h; % function handle for measurement

for k = 1:maxk
   fprintf('k = %d\n',k);
   Y{k} = interval(system.z(:,k) + system.Wsim.bounds)';
   psv{k} = PS_vsivia(Y{k},U0);
   psv{k}.epsilon = '0.01%';
   [Sz{k},Ez{k},Nz{k}] = vsivia(psv{k});
%    P = propagate(S{k}); % propagate solution boxes
%    U0 = P.box() + padding; % set new start interval for vsivia
%    intersect = true; % activate intersection mode for next iteration
end