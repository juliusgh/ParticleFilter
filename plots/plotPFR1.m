figure(1);
plotIntervals(S{21},'g',1);
plotIntervals(E{21},'k',1);
hold on
PFR.plotParticles(20);
PFR.plotX(20);
xlim([-1 2]); ylim([1 4]);
hold off
fig2Eps4Latex