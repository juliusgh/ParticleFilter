system = System();
alg = PFAlgorithms.WeightOptimisation;
options.iterations = 16;
options.debug = false;
options.gurobi = false;
sampleList = [];
timeList = [];
for i = 1:30
    options.samples = i*100;
    time = testPF(system,alg,options);
    sampleList(end+1) = options.samples;
    timeList(end+1) = time;
    fprintf('samples: %d, time: %f\n',options.samples,time);
end

function time = testPF(system,alg,options)
    rounds = 10;
    times = NaN(1,rounds);
    for i=1:rounds
        PF = ParticleFilter(system,alg,options);
        PF.run();
        times(i) = PF.time;
    end
    time = mean(times);
end