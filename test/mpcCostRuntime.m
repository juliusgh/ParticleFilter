Ns = [1;10;100;1000;10000;50000;100000];
P = PF.particlesAt(5);
I = [S{4} E{4}];
for j=1:2
    I = [I I];
end
for j=1:1
    P = [P P];
end
u = zeros(2,mpc.H);
times1 = NaN(1,numel(Ns));
times2 = NaN(1,numel(Ns));
times3 = NaN(1,numel(Ns));
times4 = NaN(1,numel(Ns));
xref = zeros(mpc.system.dim,mpc.H+1);
uref = mpc.system.uref(xref);

for i = 1:numel(Ns)
    N = Ns(i);
    PN = P(:,1:N);
    IN = I(1:N,:);
    
    tic1 = tic();
    PNt = mpc.system.predict(PN,u);
    times1(i) = toc(tic1);
    
    tic2 = tic();
    INt = mpc.system.predictIntv(IN,u);
    times2(i) = toc(tic2);
    
    tic3 = tic();
    c = mpc.Jv(PNt,u,xref,uref);
    PNc = max(c);
    times3(i) = toc(tic3);
    
    tic4 = tic();
    c = mpc.JIntv(INt,u,xref,uref);
    INc = max(c.upper(:));
    times4(i) = toc(tic4);
    disp(['N=',num2str(N),' t1=',num2str(times1(i)),' t2=',num2str(times2(i)),' t3=',num2str(times3(i)),' t4=',num2str(times4(i))]);
end