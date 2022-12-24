tic

rho = 2.7e-6;
F = 685000; %N
w1 = .5;
w2 = .5;

%V_hb = @(x) pi*(x(1)^2 - x(2)^2)*x(7);
%V_ht = @(x) pi*(x(3)^2 - x(4)^2)*x(8);
%V_fork = @(x) pi*(x(5)^2 - x(6)^2)*(2*x(9) + x(10) + 2*x(11)/sind(x(12)/2));
%Mass = (V_hb + V_ht + V_fork)*rho ;
Mass = @(x)(pi*(x(1)^2 - x(2)^2)*x(7) + pi*(x(3)^2 - x(4)^2)*x(8) + ...
    pi*(x(5)^2 - x(6)^2)*(x(9) + 2*x(10) + 2*x(11)/sind(x(12)/2)))*rho;
Cycle  = @(x) -(14479/((F/(pi*(x(3)^2 - x(4)^2) + 2*pi*(x(5)^2 - x(6)^2))) ...
    - 96.5))^2;

f = @(x) w1*x(14) + w2*x(15);
lb = [30, 27, 33, 30, 35, 33, 680, 350, 100, 335, 40, 71.5, 0, 0, 0, 0];
ub = [33, 30, 35, 32, 40, 38, 730, 500, 300, 430, 55, 74.5, 0, 100, 10000000, 0];
A = [-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
      0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
      0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
b = [-2,-2,-2];
Aeq = [];
beq = [];

x0 = [31, 29, 34, 29, 39, 37, 729, 499, 299, 429, 54, 74.4, 0, 0, 0, 0];


history= histclass;
outf= @(x,optimValues,state)outfun3(x,optimValues,state,history);

nonlcon = @nonlincon;
options = optimset('OutputFcn',outf,'Algorithm','interior-point','PlotFcn', ...
    {@optimplotx,@optimplotfval,@optimplotfirstorderopt});

[x,output,history,lambda] = fmincon(f,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)
Mass(x)
-Cycle(x)


toc
