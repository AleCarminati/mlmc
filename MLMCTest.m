%%
% Reproduce the top left and top right plots in Figure 2 in the article.
clc
rng(248963849)
mlmc = MLMC(1,16,800,1,0.3,1);
mlmc.plot_var_mean(5,10000);

%%
% Reproduce the top left and top right plots in Figure 3 in the article.
clc
rng(80392666)
mlmc = MLMC(2,8,1400,1,0.3,1);
mlmc.plot_var_mean(4,1000);

%%
% Run the MLMC to reproduce the bottom plots as in the bottom part of 
% Figure 2.  
clc
% Vector of epsilon values for the 1D case.
epsilon_vec = [0.01, 0.0075, 0.005, 0.0025];
% Vector of epsilon values for the 2D case.
% epsilon_vec = [0.01, 0.0075, 0.005, 0.0025];
count = 1;
costs = zeros('like',epsilon_vec);
n_samples = zeros(size(epsilon_vec,2), 5);
for epsilon=epsilon_vec
	 mlmc = MLMC(1,16,800,1,0.3,1);
	% mlmc = MLMC(2,8,1400,1,0.3,1);
	mlmc = mlmc.run_epsilon_fixed(1.75,2,epsilon);
	costs(count) = mlmc.computeCost();
	n_levels = size(mlmc.levels,1);
	for idx=1:n_levels
		n_samples(count, idx)=mlmc.levels(idx).N_l;
	end
	count = count+1;
end

%%
% Plot the number of samples in each level.
clc

for count=1:size(epsilon_vec,2)
	semilogy(0:4,n_samples(count,:), "-o", ...
		'DisplayName',"$\varepsilon=$"+num2str(epsilon_vec(count)))
	hold on
end

legend("Interpreter","latex")

xticks(0:(n_levels-1))
xlabel("Level")
ylabel("$N_\ell$","Interpreter","latex")
hold off

%%
% Plot the cost for each MLMC run.
loglog(epsilon_vec, costs.*epsilon_vec.^2, "-o")
ylim([10,100])
xlabel("Accuracy $\varepsilon$","Interpreter","latex")
ylabel("$\varepsilon^2$ cost","Interpreter","latex")