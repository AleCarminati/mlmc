% This script tests the class FVSolver, comparing graphically its solutions 
% with the correct ones.

%%
clc 
d = 1;
m = 256;
k = @(x) sqrt(x)+1;
solver = FVSolver(d,@(pointsSet,m) kWrapper(1,k,pointsSet,m),m);

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
solver = FVSolver(d,@(pointsSet,m) kWrapper(1,k,pointsSet,m),m);

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
solver = FVSolver(d,@(pointsSet,m) kWrapper(1,k,pointsSet,m),m);

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
solver = FVSolver(d,@(pointsSet,m) kWrapper(1,k,pointsSet,m),m);

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
solver = FVSolver(d,@(pointsSet,m) kWrapper(1,k,pointsSet,m),m);

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
solver = FVSolver(d,@(pointsSet,m) kWrapper(1,k,pointsSet,m),m);

x_plot = linspace(0,1,1024);
y_plot = solver.getSolutionValue(x_plot);
f_true = @(x) (-x-exp(x)+exp(1)+1)/exp(1);
y_true = f_true(x_plot);
plot(x_plot, y_plot,x_plot, y_true)

%%
clc 
d = 2;
m = 16;
k = @(x,y) sinh(x)+exp(-15.*(1-y));
solver = FVSolver(d,@(pointsSet,m) kWrapper(2,k,pointsSet,m),m);

[x_plot,y_plot] = meshgrid(linspace(0,1,1024));
z_plot = solver.getSolutionValue(x_plot,y_plot);
surf(x_plot,y_plot,z_plot)

%%
clc 
d = 2;
m = 32;
k = @(x,y) sinh(0.5*(x+1)./(y+1));
solver = FVSolver(d,@(pointsSet,m) kWrapper(2,k,pointsSet,m),m);

[x_plot,y_plot] = meshgrid(linspace(0,1,1024));
z_plot = solver.getSolutionValue(x_plot,y_plot);
surf(x_plot,y_plot,z_plot)

%%
% This function creates a wrapper for the functions in input, so that they
% are compatible with FVSolver.
function value = kWrapper(d, f, varargin)
	pointsSet = varargin{1};
	m = varargin{2};	
	if d==1
		switch pointsSet
			case "cen"
				xpoints = linspace(1/m,1-1/m,m-1)';
			case "ext"
				xpoints = [0,1]';
		end
		value = f(xpoints);
	else
		switch pointsSet
			case "cpx"
				[xpoints, ypoints] =...
					meshgrid(linspace(1/(2*m),1-1/(2*m),m), ...
					linspace(1/(m),1-1/(m),m-1));
				xpoints = reshape(xpoints,[],1);
				ypoints = reshape(ypoints,[],1);
			case "cpy"
				[xpoints, ypoints] =...
					meshgrid(linspace(1/(m),1-1/(m),m-1), ...
					linspace(1/(2*m),1-1/(2*m),m));
				xpoints = reshape(xpoints,[],1);
				ypoints = reshape(ypoints,[],1);
			case "bl"
				xpoints =repelem(0,m)';
				ypoints = linspace(1/(2*m),1-1/(2*m),m)';
			case "br"
				xpoints =repelem(1,m)';
				ypoints = linspace(1/(2*m),1-1/(2*m),m)';
		end
		value = f(xpoints,ypoints);
	end
end