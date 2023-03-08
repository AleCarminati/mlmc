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

			variances = zeros(n_levels,1);
			means = zeros(n_levels,1);
			for i=1:n_levels
				fprintf("Level %d\n", i);
				level = Level(100, obj.d, obj.m_0*2^i, false, obj.rfs);
				means(i) = level.Y_l;
				variances(i) = var(level.Y_vec);
			end

			figure(1)
			plot(1:n_levels,log2(abs(means)), "*-")
			title("Mean of $Q_l-Q_{l-1}$", "Interpreter","latex")
			xticks(1:n_levels)
			xlabel("Level")
			ylabel("$\log_2|$mean$|$","Interpreter","latex")
			figure(2)
			plot(1:n_levels, log2(variances), "*-")
			title("Variance of $Q_l-Q_{l-1}$", "Interpreter","latex")
			xticks(1:n_levels)
			xlabel("Level")
			ylabel("$\log_2$ variance","Interpreter","latex")
		end

		function obj = run_epsilon_fixed(obj, alpha, eps)
			% This function runs the MLMC algorithm for a fixed error tolerance
			% eps and a fixed convergence rate alpha.

			N_0 = 100;
			obj.levels = Level(N_0,obj.d, obj.m_0, true,obj.rfs);
			n_levels = 1;
			convergence = false;
			while not(convergence)
				disp(n_levels)
				% Compute optimal number of samples for each level and evaluate 
				% extra samples
				for idx=1:n_levels
					var_level = var(obj.levels(idx).Y_vec);
					cost_level = (obj.m_0*2^(idx-1))^obj.gamma;
					obj.levels(idx) = obj.levels(idx).updateNumSamples( ...
						ceil(sqrt(var_level/cost_level)));
				end
	
				% Convergence test
				if n_levels>1
					max_mean_diff = 0;
					for idx=n_levels:max(n_levels-2,2)
						if(abs(obj.levels(idx).Y_l)>max_mean_diff)
							max_mean_diff = abs(obj.levels(idx).Y_l);
						end
					end
					if (max_mean_diff/(2^alpha-1)<eps/sqrt(2))
						convergence = true;
					end
				end

				if not(convergence)
					obj.levels = [obj.levels; 
							Level(N_0,obj.d, obj.m_0*2^n_levels, true,obj.rfs)];
						n_levels = n_levels + 1;
				end
			end
		end
	end
end