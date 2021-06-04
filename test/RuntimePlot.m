%load('data/RuntimeTest.mat')
figure(1)
plot(sampleList_N,timeList_N,'*-',sampleList_R,timeList_R,'*-',sampleList_W,timeList_W,'*-',sampleList_WG,timeList_WG,'*-')
ylim([0 10])
xlim([1000 10000])
xlabel('Anzahl an Partikeln')
ylabel('Laufzeit in s')
%title('Algorithmen-Vergleich: Laufzeit bzgl. Partikelzahl')
legend('Naive PF (Alg. 3)','Red. Remov. PF (Alg. 4)','Weight Opt. PF (Alg. 5)','Weight Opt. PF (Alg. 5, GUROBI)','Location','north')
itm_formatfig(2)

% figure(2)
% plot(sampleList_W,timeList_W,'*-',sampleList_WG,timeList_WG,'*-')
% C = colororder(gca);
% C = C(3:end,:);
% colororder(gca,C);
% ylim([0 170])
% xlim([1000 10000])
% xlabel('Anzahl an Partikeln')
% ylabel('Laufzeit in s')
% title('Algorithmen-Vergleich: Laufzeit bzgl. Partikelzahl')
% legend('Weight Opt. PF','Weight Opt. PF (GUROBI)','Location','southeast')
