figure(1);
k = 20;
plotIntervals(PFN.intvS{k+1},'g',1);
%plotIntervals(PFW.intvE{k+1},'g',1);
hold on
PFW.plotParticles(k);
%PFCC.plotConvhull(k);
PFW.plotX(k);
xlim([-1 2]);
ylim([1 4]);
xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
hold off
itm_formatfig(4,'FigSize',[8;8])