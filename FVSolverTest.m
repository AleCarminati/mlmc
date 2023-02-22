% This script tests the class FVSolver, comparing graphically its solutions 
% with the correct ones.

%%
clc 
d = 1;
m = 1024;
k = @(x) sqrt(x)+1;
solver = FVSolver(d,k,m);

x_plot = linspace(0,1,1024);
y_plot = solver.getSolutionValue(x_plot);
f_true = @(x) (sqrt(x)-log(sqrt(x)+1)-1+log(2))/(log(2)-1);
y_true = f_true(x_plot);
plot(x_plot, y_plot,x_plot, y_true)

%%
clc 
d = 1;
m = 256;
k = @(x) cosh(x);
solver = FVSolver(d,k,m);

x_plot = linspace(0,1,1024);
y_plot = solver.getSolutionValue(x_plot);
f_true = @(x) 1-(atan(tanh(x/2)))/(atan(tanh(1/2)));
y_true = f_true(x_plot);
plot(x_plot, y_plot,x_plot, y_true)

%%
clc 
d = 1;
m = 256;
k = @(x) exp(x);
solver = FVSolver(d,k,m);

x_plot = linspace(0,1,1024);
y_plot = solver.getSolutionValue(x_plot);
f_true = @(x) exp(-x).*(exp(1)-exp(x))/(exp(1)-1);
y_true = f_true(x_plot);
plot(x_plot, y_plot,x_plot, y_true)

%%
clc 
d =1;
m = 256;
k = @(x) 1+x.*0;
solver = FVSolver(d,k,m);

x_plot = linspace(0,1,1024);
y_plot = solver.getSolutionValue(x_plot);
f_true = @(x) 1-x;
y_true = f_true(x_plot);
plot(x_plot, y_plot,x_plot,y_true)

%%
clc 
d = 1;
m = 256;
k = @(x) x.^2+x+2;
solver = FVSolver(d,k,m);

x_plot = linspace(0,1,1024);
y_plot = solver.getSolutionValue(x_plot);
f_true = @(x) (atan((2*x+1)/sqrt(7))-atan(3/sqrt(7)))./ ... 
	(atan(1/sqrt(7))-atan(3/sqrt(7)));
y_true = f_true(x_plot);
plot(x_plot, y_plot,x_plot, y_true)

%%
clc 
d = 1;
m = 256;
k = @(x) 1./(1+exp(x));
solver = FVSolver(d,k,m);

x_plot = linspace(0,1,1024);
y_plot = solver.getSolutionValue(x_plot);
f_true = @(x) (-x-exp(x)+exp(1)+1)/exp(1);
y_true = f_true(x_plot);
plot(x_plot, y_plot,x_plot, y_true)

%%
clc 
d = 1;
m = 256;
k = @(x) x.^(-2);
solver = FVSolver(d,k,m);

x_plot = linspace(0,1,1024);
y_plot = solver.getSolutionValue(x_plot);
f_true = @(x) 1-x.^3;
y_true = f_true(x_plot);
plot(x_plot, y_plot,x_plot, y_true)