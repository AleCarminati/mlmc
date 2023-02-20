d =1;
m = 16;
k = @(x)sin(x);
solver = FVSolver(d,k,m);
solver.getSolutionValue(0.35);

%%
x_plot = linspace(0,1);
