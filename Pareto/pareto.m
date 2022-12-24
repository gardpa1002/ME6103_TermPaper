rho = 2.7e-6;
F = 685000; %N
tic
%V_hb = @(x) pi*(x(1)^2 - x(2)^2)*x(7);
%V_ht = @(x) pi*(x(3)^2 - x(4)^2)*x(8);
%V_fork = @(x) pi*(x(5)^2 - x(6)^2)*(2*x(9) + x(10) + 2*x(11)/sind(x(12)/2));
%Mass = (V_hb + V_ht + V_fork)*rho ;
Mass = @(x)(pi*(x(1)^2 - x(2)^2)*x(7) + pi*(x(3)^2 - x(4)^2)*x(8) + pi*(x(5)^2 - x(6)^2)*(x(9) + 2*x(10) + 2*x(11)/sind(x(12)/2)))*rho;
Cycle  = @(x) -(14479/((F/(pi*(x(3)^2 - x(4)^2) + 2*pi*(x(5)^2 - x(6)^2))) - 96.5))^2;


lb = [30, 27, 33, 30, 35, 33, 680, 350, 100, 335, 40, 71.5];
ub = [33, 30, 35, 32, 40, 38, 730, 500, 300, 430, 55, 74.5];
A = [-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
      0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0;
      0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0];
b = [-2,-2,-2];
Aeq = [];
beq = [];



fun = @(x) [(pi*(x(1)^2 - x(2)^2)*x(7) + pi*(x(3)^2 - x(4)^2)*x(8) + pi*(x(5)^2 - x(6)^2)*(x(9) + 2*x(10) + 2*x(11)/sind(x(12)/2)))*rho;
    -(14479/((F/(pi*(x(3)^2 - x(4)^2) + 2*pi*(x(5)^2 - x(6)^2))) - 96.5))^2];


nonlcon = @nonlinconpareto;
opts_ps = optimoptions('paretosearch','Display','off','PlotFcn','psplotparetof');
rng default % For reproducibility

[x,fval,exitflag,output,residuals] = paretosearch(fun,12,A,b,Aeq,beq,lb,ub,[],opts_ps)
toc