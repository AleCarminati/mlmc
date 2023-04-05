%%
% Test that the computed eigenvalues are as the ones represented in Figure
% 1 in the article.
clc
r = RandomFieldSampler(1000,1,0.01,1);
loglog(1:length(r.eigenvalues),r.eigenvalues)
hold on
r = RandomFieldSampler(1000,1,0.1,1);
loglog(1:length(r.eigenvalues),r.eigenvalues,"--")
hold on
r = RandomFieldSampler(1000,1,1,1);
loglog(1:length(r.eigenvalues),r.eigenvalues,"-.")
legend("\lambda=0.01","\lambda=0.1","\lambda=1")
ylim([1e-8 1])

%% Test if w_n are the solution of the transcendental equation
clc
m_kl = 1400;
lambda = 0.3;
r = RandomFieldSampler(m_kl,1,lambda,1);
sol = abs(tan(r.w)-(2*lambda.*r.w)./(lambda^2.*r.w.^2-1));
max(sol)

%% Graphically test if w_n are the solution of the transcendental equation
clc
m_kl = 1400;
lambda = 0.3;
r = RandomFieldSampler(m_kl,1,lambda,1);
x_plot = linspace(0,max(r.w),10000);
y_tan = tan(x_plot);
y_frac = 2*lambda*x_plot./(lambda^2.*x_plot.^2-1);
% These instructions remove the vertical lines in the plot of the function
y_tan(diff([0,y_tan])<0)=NaN;
y_frac(diff([0,y_frac])>0)=NaN;
plot(x_plot,y_tan)
hold on
plot(x_plot,y_frac)
plot(r.w,tan(r.w),".","color","green", "MarkerSize",15)
ylim([-10,10])

%% Test if A_n are such that the L-2 norm of b_n is 1.
clc
n_int = 10000;
m_kl = 1400;
r = RandomFieldSampler(m_kl,1,0.3,1);
x_int = linspace(0,1,n_int);
vec_int = zeros(m_kl,1);
for idx=1:n_int
	vec_int = vec_int + r.b_n(x_int(idx)).^2/n_int;
end

max(abs(vec_int-1))

%% Create and plot a 1D random field
clc
m_kl = 800;
sigma2 = 1;
lambda = 0.3;
n_points_plot = 1000;
r = RandomFieldSampler(m_kl,sigma2,lambda,1);
x_plot = linspace(0,1,n_points_plot)';
y_plot = r.computeRandomFieldValue(x_plot);
plot(x_plot,y_plot)

%% Create and plot a 2D random field
clc
m_kl = 1400;
sigma2 = 1;
lambda = 0.3;
n_points_plot = 100;
r = RandomFieldSampler(m_kl,sigma2,lambda,2);
[x_plot, y_plot] = meshgrid(linspace(0,1,n_points_plot));
x_plot_vec = reshape(x_plot,[],1);
y_plot_vec = reshape(y_plot,[],1);
z_plot_vec = r.computeRandomFieldValue(x_plot_vec,y_plot_vec);
z_plot = reshape(z_plot_vec,n_points_plot,n_points_plot);
surf(x_plot,y_plot,z_plot)