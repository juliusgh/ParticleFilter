fig = 1;
k = 20;
figure(fig);
plotIntervals(S{k+1},'g',fig);
plotIntervals(E{k+1},'g',fig);
hold on
PFN.plotX(k);
hold off
xlabel('$x_1$','Interpreter','latex');
ylabel('$x_2$','Interpreter','latex');
xlim([-1.5 2]); ylim([1 4.5]);
itm_formatfig(1);