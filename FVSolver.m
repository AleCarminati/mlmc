classdef FVSolver	
	properties
		d % Number of dimensions.
		m % Mesh size.
		solutionPoints % Matrix (or vector) containing the value of the 
			% solution in each element of the grid.   
	end
	
	methods
		function obj = FVSolver(d,k,m)
			obj.d = d;
			obj.m = m;
			if d==1
				if nargin(k)~=1
					error("Invalid function k, it requires %d inputs while it"+...
						"should require 1 input", nargin(k));
				end
				obj.solutionPoints = zeros(m);
			elseif d==2
				if nargin(k)~=2
					error("Invalid function k, it requires %d inputs while it"+...
						"should require 2 inputs", nargin(k));
				end
				obj.solutionPoints = zeros(m^2);
			else
				error("Invalid input: d=%d, but d must be equal to 1 or 2.", d);
			end
			
			if d==1
				k0 = k(0);
				k1 = k(1);
				centralPointsGrid = linspace(1/m,1-1/m,m);
				kCentralPointsGrid = k(centralPointsGrid);
				A = zeros(m,m);
				b = zeros(m,1);
				b(1) = -2*k0;
				harmMeans = 2./(1./kCentralPointsGrid(1:m-1)+1./ ...
					kCentralPointsGrid(2:m));
				A = A+diag(harmMeans,1)+diag(harmMeans,-1);
				A = A-diag([2*k0+harmMeans(1),harmMeans(1:m-2)+harmMeans(2:m-1),...
					2*k1+harmMeans(m-1)]);
				obj.solutionPoints = A\b;
			else
				%
			end
		end

		function value = getSolutionValue(obj,x)
			if obj.d==1
				if not(all(x>=0)) || not(all(x<=1))
					error("Invalid input: x must be in the range [0,1].");
				end 
				value = obj.solutionPoints(max(1,ceil(x*obj.m)));
			end
		end
	end
end

