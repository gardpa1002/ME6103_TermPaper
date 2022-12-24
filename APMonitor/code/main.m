addpath('../apm');

% Select server
server = 'http://byu.apmonitor.com';

% Application name
app = 'trial';

% Clear previous application
apm(server,app,'clear all');

% Load model file
apm_load(server,app,'multi_obj.apm');

% Option to select solver (1=APOPT, 2=BPOPT, 3=IPOPT)
apm_option(server,app,'nlc.solver',1);

% Solve on APM server
apm(server,app,'solve')

% Retrieve results
disp('Results');
results = apm_sol(server,app)
values = cell2mat(results(2:end,:));

% Display Results in Web Viewer 
url = apm_web_var(server,app);
	
%x = results.values(7:18)
x = results.values(11:26);
rho = 2.7e-6;
F = 685000; %N

Mass = @(x)(pi*(x(1)^2 - x(2)^2)*x(7) + pi*(x(3)^2 - x(4)^2)*x(8) + pi*(x(5)^2 - x(6)^2)*(x(9) + 2*x(10) + 2*x(11)/sin(x(12)/2)))*rho;
Cycle  = @(x) (14479/((F/(pi*(x(3)^2 - x(4)^2) + 2*pi*(x(5)^2 - x(6)^2))) - 96.5))^2;

Mass(x)
Cycle(x)
