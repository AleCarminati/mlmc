%%
% Reproduce the top left and top right plots in Figure 2 in the article.
clc
rng(828782066)
mlmc = MLMC(1,16,800,1,0.3,1);
mlmc.plot_var_mean(5,100000);

%%
% Reproduce the top left and top right plots in Figure 3 in the article.
clc
rng(80392666)
mlmc = MLMC(2,8,1400,1,0.3,1);
mlmc.plot_var_mean(4,1000);

