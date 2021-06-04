figure(1);
plotX = subplot(2,1,1);
plotU = subplot(2,1,2);

figX = figure(2);
copyobj(plotX,figX);
xlabel('Zeitschritt $k$','Interpreter','latex')
ylabel('\textbf{\emph x}$_k$ in rad','Interpreter','latex')
legend('$x_1$','$x_2$','Interpreter','latex','NumColumns',2)
title('')
itm_formatfig(2,'LegendLocation','South','FigSize',[14;4])

figU = figure(3);
copyobj(plotU,figU);
xlabel('Zeitschritt $k$','Interpreter','latex')
ylabel('\textbf{\emph u}$_k$ in Nm','Interpreter','latex')
legend('$u_1$','$u_2$','Interpreter','latex','NumColumns',2)
title('')
itm_formatfig(2,'LegendLocation','South','FigSize',[14;4])