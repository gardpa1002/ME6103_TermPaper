function [c,ceq] = nonlinconsingle(x)
E = 55;
Mass = 100;


I = (pi/2)*(x(3)^4 - x(4)^4);
F1 = Mass*9.81*sind(x(12));

c = x(8) + x(10) + x(9) - 2*sqrt((pi*E*I)/F1);
ceq = [];

end