% parameters
tmax = 10; % timesteps
system = ManipulatorSystemDirect();

% initialize MPC
mpc0 = MPC(system,0); % 0: no filter, use measurments
mpc1 = MPC(system,1,50e3); % ParticleFilter with 50e3 particles
mpc1.simulate = false;
mpc2 = MPC(system,2); % 2: no filter, use exact system states
mpc2.simulate = false;

% run MPC
disp('run mpc0...');
tic();
mpc0.runUntil(tmax);
T0 = toc();
disp('finished mpc0!');
mpc1.v = mpc0.v;
mpc1.w = mpc0.w;
disp('run mpc1...');
tic();
mpc1.runUntil(tmax);
T1 = toc();
disp('finished mpc1!');
mpc2.v = mpc0.v;
mpc2.w = mpc0.w;
disp('run mpc2...')
tic();
mpc2.runUntil(tmax);
T2 = toc();
disp('finished mpc2!');

% evaluate MPC
J0 = mpc0.realJ();
J1 = mpc1.realJ();
J2 = mpc2.realJ();
disp(['mpc0: T = ',num2str(T0),', J = ',num2str(J0)]);
disp(['mpc1: T = ',num2str(T1),', J = ',num2str(J1)]);
disp(['mpc2: T = ',num2str(T2),', J = ',num2str(J2)]);