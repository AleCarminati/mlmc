%%
% Test that the computed eigenvalues are as the ones represented in Figure
% 1 in the article.
clc
r = RandomFieldSampler(1000,1,0.01,1);
loglog(1:length(r.eigenvalues),r.eigenvalues)
hold on
r = RandomFieldSampler(1000,1,0.1,1);
loglog(1:length(r.eigenvalues),r.eigenvalues)
hold on
r = RandomFieldSampler(1000,1,1,1);
loglog(1:length(r.eigenvalues),r.eigenvalues)
legend("0.01","0.1","1")

%%
clc
r = RandomFieldSampler(1000,1,1,2);
[x_plot, y_plot] = meshgrid(linspace(0,1,100));
x_plot_vec = reshape(x_plot,[],1);
y_plot_vec = reshape(y_plot,[],1);
z_plot_vec = r.computeRandomFieldValue(x_plot_vec,y_plot_vec);
z_plot = reshape(z_plot_vec,100,100);
surf(x_plot,y_plot,z_plot)
