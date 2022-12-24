function [c,ceq] = nonlinconpareto(x_ps1)
E = 68;
Mass = 100;
x = x_ps1;

I = (pi/2)*(x(3)^4 - x(4)^4);
F1 = Mass*9.81*sind(x(12));

c(1) = x(8) + x(10) + x(9) - 2*sqrt((pi*E*I)/F1);
%c(2) = x(8) + x(10) + x(9) - 2*sqrt((pi*E*I)/F1);
ceq = [];

end