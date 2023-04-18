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
		d % Number of dimensions
		m % The value of the biggest mesh size used in the level. If the level 
			% is Q_{M_L}-Q_{M_{L-1}}, m is equal to M_L.
		isLevel0 % Represents if the class represents the first level.
		rfs % Object of class RandomFieldSampler to sample from the 
			% random field.
	end
	
	methods
		function obj = Level(N_l, d, m, isLevel0, randFieldSampl)
			obj.N_l = 0;
			obj.d = d;
			obj.m = m;
			obj.isLevel0 = isLevel0;
			obj.rfs = randFieldSampl;
			obj = obj.updateNumSamples(N_l);
		end
		
		function obj = updateNumSamples(obj,N_l)
			% Updates the number of samples for the level. If the new number of
			% samples is greater than the old one, it generates new samples. 
			% Otherwise, it does nothing.

			if N_l>obj.N_l
				% Increase the dimensions of the vector of samples.
				temp_vector = zeros(N_l,1);
				temp_vector(1:obj.N_l) = obj.Y_vec;
				obj.Y_vec = temp_vector;

				diffN = N_l-obj.N_l;
				for i=1:diffN
					sampled_value = obj.getNewSample();
					obj.Y_vec(obj.N_l + i) = sampled_value; 
				end
				obj.N_l = N_l;
			end
		end
	end
	
	methods(Access=private)
		function value = getNewSample(obj)
			% Generates a new sample for the level and computes the value of the
			% QoI using an FVSolver object. If the level is not the first level,
			% it also subtracts the QoI computed from the previous level.

			obj.rfs = obj.rfs.updateRandom();
			% Workaround to pass the method as a function handle.
			if obj.d==1
				f = @(x) obj.rfs.computeRandomFieldValue(x);
			else
				f = @(x,y) obj.rfs.computeRandomFieldValue(x,y);
			end
			solver = FVSolver(obj.d, f,obj.m);
			value = obj.getQoI(solver);
			if not(obj.isLevel0)
				solver = FVSolver(obj.d, f,obj.m/2);
				value = value-obj.getQoI(solver);
			end
		end

		function value = getQoI(obj,solver)
			% Computes the value of the QoI for the given solver object. 

			if obj.d == 1 
				derivative = ... 
					(solver.solutionPoints(end)-solver.solutionPoints(end-1))* ...
					solver.m;
				value = - obj.rfs.computeRandomFieldValue( ...
					(solver.m-1)/solver.m)*derivative;

			else
				idx =(1:solver.m)';
				derivative = ...
 						(0-solver.solutionPoints(solver.m*idx))*solver.m*2;
				value = - 1/solver.m* ...
 						obj.rfs.computeRandomFieldValue( ...
						repelem(1,solver.m)',idx/solver.m-1/(2*solver.m))'*derivative;
			end
		end
	end
end

