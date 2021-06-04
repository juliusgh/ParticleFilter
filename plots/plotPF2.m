figure(1);
P = PFW.particlesAt(20);
histogram(P(1,:),'FaceColor','b')
xlabel('$x_1$','Interpreter','latex')
ylabel('$\mathrm{Anzahl~Partikel}$','Interpreter','latex')
itm_formatfig(1)