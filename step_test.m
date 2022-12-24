%----------------------------------------------------------------------------------------------------
% An implementation of Steepest descent method for optimization problems
% https://github.com/EiriniMits/Optimization-Techniques/blob/master/Steepest_Descent.m
%----------------------------------------------------------------------------------------------------

tic
clear all;
close all;
format short

N=12;
xi=sym(zeros(1,N));
for i=1:N 
    syms("x"+i);
    xi(1,i) = ("x"+i);
end

yi=sym(zeros(1,N));
for i=1:N 
    syms("y"+i);
    yi(1,i) = ("y"+i);
end
rho = 2.7e-6;

%V_hb = @(x) pi*(x(1)^2 - x(2)^2)*x(7);
%V_ht = @(x) pi*(x(3)^2 - x(4)^2)*x(8);
%V_fork = @(x) pi*(x(5)^2 - x(6)^2)*(2*x(9) + x(10) + 2*x(11)/sind(x(12)/2));
%Mass = (V_hb + V_ht + V_fork)*rho ;
Mass = @(x)(pi*(x(1)^2 - x(2)^2)*x(7) + pi*(x(3)^2 - x(4)^2)*x(8) + pi*(x(5)^2 - x(6)^2)*(x(9) + 2*x(10) + 2*x(11)/sind(x(12)/2)))*rho;
Cycle  = @(x) -(14479/((F/(pi*(x(3)^2 - x(4)^2) + 2*pi*(x(5)^2 - x(6)^2))) - 96.5))^2;





x_old = zeros(1,12);
x = [31, 29, 34, 29, 39, 37, 729, 499, 299, 429, 54, 74.4];% Initial Guess
e = 10^(-5); % Convergence Criteria
k = 1; % Iteration Counter
r = 10000;
c
new = 1;
old = 0;
E = 55;
Mass = 100;
rho = 2.7e-6;
% F = 685000; %N
% I = (pi/2)*(x3^4 - x4^4);
% F1 = Mass*9.81*sind(x12);
% const = x8 + x10 + x9 - 2*sqrt((pi*E*I)/F1);
% g = -1/const;
% obj = (pi*(x1^2 - x2^2)*x7 + pi*(x3^2 - x4^2)*x8 + pi*(x5^2 - x6^2)*(x9 + 2*x10 + 2*x11/sind(x12/2)))*rho;
% f = obj + r*g;

%%%%%%%%%%%%%%%%%
y1 = 30 + (33-30)*sin(x1).^2;
y2 = 27 + (30-27)*sin(x2).^2;
y3 = 33 + (35-33)*sin(x3).^2;
y4 = 30 + (32-30)*sin(x4).^2;
y5 = 35 + (40-35)*sin(x5).^2;
y6 = 33 + (38-33)*sin(x6).^2;
y7 = 680 + (730-680)*sin(x7).^2;
y8 = 350 + (500-350)*sin(x8).^2;
y9 = 100 + (300-100)*sin(x9).^2;
y10 = 335 + (430-335)*sin(x10).^2;
y11 = 40 + (55-40)*sin(x11).^2;
y12 = 71.5 + (74.5-71.5)*sin(x12).^2;

I = (pi/2)*(y3^4 - y4^4);
F1 = Mass*9.81*sind(y12);
const = y8 + y10 + y9 - 2*sqrt((pi*E*I)/F1);
const1 = -y1 + y2 +2;
const2 = -y3 + y4 +2;
const3 = -y5 + y6 +2;
g = -((1/const) +  (1/const1)+(1/const2)+(1/const3));
obj = (pi*(y1^2 - y2^2)*y7 + pi*(y3^2 - y4^2)*y8 + pi*(y5^2 - y6^2)*(y9 + 2*y10 + 2*y11/sind(y12/2)))*rho;
f = obj + r*g;
%%%%%%%%%%%%%%%%%%


% Gradient Computation:
for i=1:N 
    df_dx(i)= diff(f, "x"+i);
end
G = subs(df_dx, xi, x);


% display table
fprintf('k: %d\t', k);
for i=1:N
   fprintf('x%d: %d\t\t\t',i, x(i)); 
end
fprintf('f(xk): %d\t\t', subs(obj,xi,x));
fprintf('||∇f(xk)||: %d\t', norm(G)); 
fprintf('\n');

%optimization routine
%while abs(x(1) - x_old(1)) > e ||  abs(x(2) - x_old(2)) > e
while abs(new - old) > e
    syms alpha; 
        
    f = obj - r*g;
    % Gradient Computation:
    for i=1:N 
        df_dx(i)= diff(f, "x"+i);
    end
    
    G = subs(df_dx, xi, x);
    Pk = -(G); % Search Direction

    x_old = x;
    for i=1:N 
      l(i)= x(i) + alpha * Pk(i);
    end
    
    for i=1:N
      x_f(i)= subs(f, xi, l); 
    end
    
    GradF = diff(x_f, alpha);
    %alpha = solve(GradF(1), alpha); % Step size
    for i=1:N 
      x(i)= x(i) + alpha * Pk(i); % Updated xk value
    end
   
    G = subs(df_dx, xi, x);
    Pk = -(G); % New Search Direction
    k = k+1;
    r = r*c;

    % Gradient Computation:
    for i=1:N 
        df_dx(i)= diff(obj, "x"+i);
    end
    f_grad = subs(df_dx, xi, x);

    for i=1:N 
        dg_dx(i)= diff(g, "x"+i);
    end
    g_grad = double(subs(dg_dx, xi, x));
    
 
    % display table
    fprintf('k: %d\t', k); 
    for i=1:N
        fprintf('x%d: %d\t',i, x(i)); 
    end
    fprintf('f(xk): %d\t', double(subs(obj,[xi], [x])));
    fprintf('\n');
    new = abs(subs(f,xi,x));
    old = abs(subs(f,xi,x_old));
    
end

% Gradient Computation:
for i=1:N 
    df_dx(i)= diff(obj, "x"+i);
end
f_grad = double(subs(df_dx, xi, x))

% Gradient Computation:
for i=1:N 
    dg_dx(i)= diff(g, "x"+i);
end
g_grad = double(subs(dg_dx, xi, x))


lamda = double(f_grad())/double(g_grad())

toc
