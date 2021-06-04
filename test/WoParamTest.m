system = System();
alg = PFAlgorithms.WeightOptimisation;
options.samples = 2000;
options.iterations = 16;
options.debug = false;
OTolList = [];
CTolList = [];
timeList = [];

% variiere OptimalityTolerance
for i = 1:1
    options.LPOptions.OptimalityTolerance = 1e-7*10^i;
    %options.LPOptions.OptimalityTolerance = 1e-2;
    %options.LPOptions.ConstraintTolerance = 1e-2/10^i;
    %options.LPOptions.Preprocess = 'basic';
    options.LPOptions.Display = 'final';
    options.LPOptions.Algorithm = 'dual-simplex';
    time = testPF(system,alg,options);
    OTolList(end+1) = options.LPOptions.OptimalityTolerance;
    %CTolList(end+1) = options.LPOptions.ConstraintTolerance;
    timeList(end+1) = time;
    fprintf('OTol: %d, time: %f\n',OTolList(end),time);
    %fprintf('CTol: %d, time: %f\n',CTolList(end),time);
end

function time = testPF(system,alg,options)
    rounds = 3;
    times = NaN(1,rounds);
    for i=1:rounds
        PF = ParticleFilter(system,alg,options);
        PF.run();
        times(i) = PF.time;
    end
    time = mean(times);
end