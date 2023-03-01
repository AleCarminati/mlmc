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

%% Create and plot a 1D random field
clc
m_kl = 800;
sigma2 = 1;
lambda = 0.01;
n_points_plot = 1000;
r = RandomFieldSampler(m_kl,sigma2,lambda,1);
x_plot = linspace(0,1,n_points_plot)';
y_plot = r.computeRandomFieldValue(x_plot);
plot(x_plot,y_plot)


%% Create and plot a 2D random field
clc
m_kl = 1400;
sigma2 = 1;
lambda = 0.01;
n_points_plot = 100;
r = RandomFieldSampler(m_kl,sigma2,lambda,2);
[x_plot, y_plot] = meshgrid(linspace(0,1,n_points_plot));
x_plot_vec = reshape(x_plot,[],1);
y_plot_vec = reshape(y_plot,[],1);
z_plot_vec = r.computeRandomFieldValue(x_plot_vec,y_plot_vec);
z_plot = reshape(z_plot_vec,n_points_plot,n_points_plot);
surf(x_plot,y_plot,z_plot)