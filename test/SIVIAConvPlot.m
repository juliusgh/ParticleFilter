% semilogx(epsilons3,volumes3,'*-')
semilogx(E,V(:,2:end),'*-','Color',[.5 .5 .5])
hold on
semilogx(E,V(:,21),'*-','Color','b')
set(gca,'XDir','reverse')
xlabel('Genauigkeit $\varepsilon$','Interpreter','latex')
ylabel('Lebesgue-Ma\ss\quad$\lambda(\mathcal{X}_{k|k})$','Interpreter','latex')
%ylim([0 2.5])
h = [plot(NaN,NaN,'*-','Color',[.5 .5 .5]);plot(NaN,NaN,'*-','Color','b')];
legend(h,'$1 \leq k \leq 19$','$k=20$','Interpreter','latex','NumColumns',2)
itm_formatfig(2)