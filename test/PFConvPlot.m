semilogx(countN,volumeN,'*-',countR,volumeR,'*-',countW,volumeW,'*-')
xlabel('Anzahl an Partikeln $N$','Interpreter','latex')
ylabel('Flächeninhalt der konvexen Hülle','Interpreter','latex')
ylim([0 1.2])
xlim([2e3 1e5])
yline(0.9458,'-.b','$\lambda(\mathcal{X}_{k|k})$','Interpreter','latex')
legend({'Naive PF (Alg. 3)','Red. Remov. PF (Alg. 4)','Weight Opt. PF (Alg. 5)'})
% C = colororder(gca);
% C = C(3:end,:);
% colororder(gca,C);

itm_formatfig(2,'LegendLocation','NorthWest')