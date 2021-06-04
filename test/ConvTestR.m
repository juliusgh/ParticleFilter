%load data/options.mat
%system = System();
alg = PFAlgorithms.RedundancyRemoval;
options.iterations = 20;
options.debug = false;

divides = 20;
samples1 = [2e3 3e3 4e3 1e4 2e4 3e4];
samples2 = [1e5 2e5 1e6 2e6 1e7 2e7];
samples = [samples1 samples2];
n1 = numel(samples1);
n2 = numel(samples2);
area = NaN(1,n1+n2);
count = NaN(1,n1+n2);
time = NaN(1,n1+n2);
for i = 1:n1
    options.samples = samples(i);
    [t, a] = testPF(system,alg,options,divides,10);
    area(i) = a;
    time(i) = t;
    count(i) = options.samples;
    fprintf('count: %d, area: %f, time: %f\n',count(i),area(i),time(i));
end
for i = (n1+1):(n1+n2)
    options.samples = samples(i);
    [t, a] = testPF(system,alg,options,divides,5);
    area(i) = a;
    time(i) = t;
    count(i) = options.samples;
    fprintf('count: %d, area: %f, time: %f\n',count(i),area(i),time(i));
end
%plot(count,area)

function [time, area] = testPF(system,alg,options,divides,rounds)
    times = NaN(1,rounds);
    areas = NaN(1,rounds);
    for i=1:rounds
        PF = ParticleFilter(system,alg,options);
        PF.run(divides);
        times(i) = PF.time;
        areas(i) = PF.area(options.iterations);
    end
    time = mean(times);
    area = mean(areas);
end