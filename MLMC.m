classdef MLMC
	
	properties
		d % Number of dimensions
		m_0 % The mesh size of the smallest level
		m_kl % Number of terms of the KL expansion
		sigma2 % Variance of exponential covariance function
		lambda % Correlation length of exponential covariance function
		gamma % Exponent of the mesh size in the cost function
		rfs % The random field sampler used in the object
		levels % Vector containing the Level classes

	end
	
	methods
		function obj = MLMC(d, m_0, m_kl, sigma2, lambda, gamma)
    	% Constructor for MLMC class
    	
			if d~=1 && d~=2
				error("Invalid input: d=%d, but d must be equal to 1 or 2.", d);
			end
			if sigma2<0
				error("Invalid input: sigma^2= %d, but the variance of the " + ...
					"exponential covariance function must be positive.");
			end

    	obj.d = d;
    	obj.m_0 = m_0;
    	obj.m_kl = m_kl;
    	obj.sigma2 = sigma2;
    	obj.lambda = lambda;
    	obj.gamma = gamma;

			obj.rfs =RandomFieldSampler(obj.m_kl, obj.sigma2, obj.lambda, obj.d);
		end

		function obj = plot_var_mean(obj, n_levels)
			% This function plots the mean and variance of each level of MLMC 
			% method, using a logarithmic scale. It takes in an arguments:
			% - n_levels: the number of levels for which to compute the mean
			%		and variance.

			N_l = 500;
			variances = zeros(n_levels,1);
			means = zeros(n_levels,1);
			for i=1:n_levels
				fprintf("Level %d\n", i);
				level = Level(N_l, obj.d, obj.m_0*2^i, false, obj.rfs);
				means(i) = mean(level.Y_vec);
				variances(i) = var(level.Y_vec);
			end

			figure(1)
			plot(1:n_levels,log2(abs(means)), "*-")
			title("Mean of $Q_l-Q_{l-1}$", "Interpreter","latex")
			xticks(1:n_levels)
			ylim([-20,0])
			yticks(-20:5:0)
			xlabel("Level")
			ylabel("$\log_2|$mean$|$","Interpreter","latex")
			figure(2)
			plot(1:n_levels, log2(variances), "*-")
			xticks(1:n_levels)
			yticks(-20:5:0)
			ylim([-20,0])
			title("Variance of $Q_l-Q_{l-1}$", "Interpreter","latex")
			xlabel("Level")
			ylabel("$\log_2$ variance","Interpreter","latex")
		end

		function obj = run_epsilon_fixed(obj, alpha, beta, eps)
			% This function runs the MLMC algorithm for a fixed error tolerance
			% eps and a fixed convergence rate alpha.

			N_0 = 20;
			obj.levels = [Level(N_0,obj.d, obj.m_0, true, obj.rfs);
				Level(N_0,obj.d, obj.m_0*2, false, obj.rfs);
				Level(N_0,obj.d, obj.m_0*2^2, false, obj.rfs)];
			n_levels = 3;
			convergence = false;
			while not(convergence)
				if n_levels ~= 3
					fprintf("\n");
				end
				fprintf("Number of levels: %d\n", n_levels)
				
				% Compute optimal number of samples for each level and evaluate 
				% extra samples.
				sum = obj.computeSumVarCost();
				for idx=1:n_levels
					var_level = var(obj.levels(idx).Y_vec);
					cost_level = (obj.m_0^obj.d*2^(idx-1))^obj.gamma;
					new_n_samples = ceil(2*eps^(-2)*sqrt(var_level/cost_level)*sum);
					fprintf("Number of samples for level %d: %d\n", ...
						idx-1,new_n_samples);
					obj.levels(idx) =obj.levels(idx).updateNumSamples(new_n_samples);
				end
	
				% Convergence test
				max_mean_diff = 0;
				for idx=n_levels:-1:(n_levels-2)
					if(abs(mean(obj.levels(idx).Y_vec))>max_mean_diff)
						max_mean_diff = abs(mean(obj.levels(idx).Y_vec));
					end
				end
				if (max_mean_diff/(2^alpha-1)<eps/sqrt(2))
					convergence = true;
				end

				% If not converged, add a new level
				if not(convergence)
					% Use the theoretical results to add to the sum also the new
					% level.
					sum = obj.computeSumVarCost();
					var_new_level = var(obj.levels(n_levels).Y_vec)/2^beta;
					cost_new_level = (obj.m_0^obj.d*2^(n_levels))^obj.gamma;
					sum = sum+sqrt(var_new_level*cost_new_level);

					n_samples = ceil(2*eps^(-2)* ...
						sqrt(var_new_level/cost_new_level)*sum); 
					obj.levels = [obj.levels; 
							Level(n_samples,obj.d, obj.m_0*2^n_levels, false,obj.rfs)];
					n_levels = n_levels + 1;
				end
			end
		end

		function cost = computeCost(obj)
			cost = 0;
			for idx=1:length(obj.levels)
				cost = cost + obj.levels(idx).N_l* ...
					obj.levels(idx).m^(obj.d*obj.gamma);
			end
		end
	end

	methods(Access=private)
		function sum = computeSumVarCost(obj)
			% Compute the sum over the levels of the square root of the product
			% of the variance and the cost, necessary to compute the optimal
			% number of samples for each level.
			sum = 0;
			for idx=1:length(obj.levels)
				var_level = var(obj.levels(idx).Y_vec);
				cost_level = (obj.m_0^obj.d*2^(idx-1))^obj.gamma;
				sum = sum+sqrt(var_level*cost_level);
			end
		end
	end
end