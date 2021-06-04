figure(1);
k = 20;
I = [S{k+1} E{k+1}];
plotIntervals(I(:,3:4),'g',1);
hold on
PFtest.plotParticles2(k);
PFtest.plotX2(k);
%xlim([-1 2]); ylim([1 4]);
xlabel('$x_3$','Interpreter','latex')
ylabel('$x_4$','Interpreter','latex')
hold off
itm_formatfig(1,'FigName','mypicture')