classdef Level
% The Level class represents a level in the multi-level Monte Carlo method
% for solving stochastic differential equations. It stores information
% about the level, such as the number of samples, the value of the
% estimator, the number of dimensions, the mesh size used in the level,
% and whether it is the first level. It also stores an object of the
% RandomFieldSampler class for sampling from the random field.

	properties
		N_l % Number of samples for that level
		Y_vec % Vector the samples for the level
		Y_l % Value of the estimator
		d % Number of dimensions
		m % The value of the biggest mesh size used in the level. If the level 
			% is Q_{M_L}-Q_{M_{L-1}}, m is equal to M_L.
		isLevel0 % Represents if the class represents the first level.
		randFieldSampl % Object of class RandomFieldSampler to sample from the 
			% random field.
	end
	
	methods
		function obj = Level(N_l, d, m, isLevel0, randFieldSampl)
			obj.N_l = 0;
			obj.Y_l = 0;
			obj.d = d;
			obj.m = m;
			obj.isLevel0 = isLevel0;
			obj.randFieldSampl = randFieldSampl;
			obj = obj.updateNumSamples(N_l);
		end
		
		function obj = updateNumSamples(obj,N_l)
			% Updates the number of samples for the level. If the new number of
			% samples is greater than the old one, it generates new samples and
			% updates the value of Y_l. Otherwise, it does nothing.

			if N_l>obj.N_l
				% Increase the dimensions of the vector of samples.
				temp_vector = zeros(N_l,1);
				temp_vector(1:obj.N_l) = obj.Y_vec;
				obj.Y_vec = temp_vector;

				partialSum = obj.N_l*obj.Y_l;
				diffN = N_l-obj.N_l;
				for i=1:diffN
					sampled_value = obj.getNewSample();
					obj.Y_vec(obj.N_l + i) = sampled_value; 
					partialSum = partialSum + sampled_value;
				end
				obj.Y_l = partialSum/N_l;
				obj.N_l = N_l;
			end
		end
	end
	
	methods(Access=private)
		function value = getNewSample(obj)
			% Generates a new sample for the level and computes the value of the
			% QoI using an FVSolver object. If the level is not the first level,
			% it also subtracts the QoI computed from the previous level.

			obj.randFieldSampl = obj.randFieldSampl.updateRandom();
			% Workaround to pass the method as a function handle.
			if obj.d==1
				f = @(x) obj.randFieldSampl.computeRandomFieldValue(x);
			else
				f = @(x,y) obj.randFieldSampl.computeRandomFieldValue(x,y);
			end
			solver = FVSolver(obj.d, f,obj.m);
			value = obj.getQoI(solver);
			if not(obj.isLevel0)
				solver = FVSolver(obj.d, f,obj.m/2);
				value = value-obj.getQoI(solver);
			end
		end

		function value = getQoI(obj,solver)
			% Computes the value of the QoI for the given solver object. For a
			% 1-dimensional problem, it numerically computes the left derivative
			% of the probability density function of the random field at x=1.
			% For a higher-dimensional problem, it numerically integrates the
			% product of the probability density function and the diffusion
			% coefficient over the domain. It uses the midpoint rule for
			% integration.

			if obj.d == 1 
				% Numerically compute the left derivative of p in 1.
				derivative = (0-solver.solutionPoints(end))*2*solver.m;
				value = - obj.randFieldSampl.computeRandomFieldValue(1)*derivative;

			else
				integral = 0;
				% For each mesh, we numerically compute the integral of p in x_1=1,
				% and we compute the integral of k using the midpoint rule for
				% integration.
				for idx=1:solver.m
					derivative = (0-solver.solutionPoints(solver.m*idx))*2*solver.m;
					integral = integral - 1/solver.m* ...
						obj.randFieldSampl.computeRandomFieldValue( ...
						1,idx/solver.m-1/(2*solver.m))*derivative;
				end
				value = integral;
			end
		end
	end
end

