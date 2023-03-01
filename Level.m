classdef Level
	properties
		N_l % Number of samples for that level
		Y_l % Value of the estimator
		d % Number of dimensions
		m % The value of the biggest mesh size used in the level. If the level 
			% is Q_{M_L}-Q_{M_{L-1}}, m is equal to M_L.
		isLevel0 % Represents if the class represents the first level.
		randFieldSampl % Object of class RandomFieldSampler to sample from the 
			% random field
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
			% If the new number of samples is greater than the old one, do
			% the update. Otherwise, do nothing.
			if N_l>obj.N_l
				partialSum = obj.N_l*obj.Y_l;
				diffN = N_l-obj.N_l;
				for i=1:diffN
					partialSum = partialSum + obj.getNewSample();
				end
				obj.Y_l = partialSum/N_l;
				obj.N_l = N_l;
			end
		end

		function value = getNewSample(obj)
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
			if obj.d == 1 
				% Numerically compute the left derivative of p in 1.
				derivative = (solver.solutionPoints(end)- ...
					solver.solutionPoints(end-1))/solver.m;
				value = - obj.randFieldSampl.computeRandomFieldValue(1)*derivative;
			else
				integral = 0;
				% For each mesh, we numerically compute the integral of p in x_1=1,
				% and we compute the integral of k using the midpoint rule for
				% integration.
				for idx=1:solver.m
					derivative = (solver.solutionPoints(solver.m*idx)- ...
						solver.solutionPoints(solver.m*idx-1))/solver.m;
					integral = integral - 1/solver.m* ...
						obj.randFieldSampl.computeRandomFieldValue( ...
						1,idx/solver.m-1/(2*solver.m))*derivative;
				end
				value = integral;
			end
		end
	end
end

