classdef FVSolver	
% This is a Matlab class called FVSolver that solves a partial differential
% equation using the Finite Volume method. It has three properties: d, 
% which is the number of dimensions of the problem; m, which is the size of
% the grid used to solve the problem; and solutionPoints, which is a matrix
% (or vector) containing the value of the solution in each element of the 
% grid.
% The class has two methods: the constructor, which takes as inputs the
% number of dimensions d, a function k, and the size of the grid m, and
% solves the problem using the Finite Volume method; and getSolutionValue,
% which takes as input a set of coordinates x and returns the value of the
% solution at those coordinates.

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
					error("Invalid function k, it requires %d inputs while it "+...
						"should require 1 input", nargin(k));
				end
				obj.solutionPoints = zeros(obj.m);
			elseif d==2
				if nargin(k)~=2
					error("Invalid function k, it requires %d inputs while it "+...
						"should require 2 inputs", nargin(k));
				end
				obj.solutionPoints = zeros(obj.m^2);
			else
				error("Invalid input: d=%d, but d must be equal to 1 or 2.", d);
			end
			
			if d==1
				A = zeros(obj.m,obj.m);
				b = zeros(obj.m,1);
				b(1) = -2*k(0);
				kValues = k(linspace(1/obj.m,1-1/obj.m,obj.m-1)');
				A = A+diag(kValues,1)+diag(kValues,-1);
				A = A-diag([2*k(0)+kValues(1); ...
					kValues(1:obj.m-2)+kValues(2:obj.m-1);...
					2*k(1)+kValues(obj.m-1)]);
				obj.solutionPoints = A\b;
			else
				midPoints_ver = reshape(k("cpy",obj.m),obj.m,obj.m-1);

				% Concatenate the midpoints of the vertical segments of the mesh 
				% grid with a column of zeros. This column is useful to put a zero 
				% every m steps, i.e. where we are in k_{i,m}.
				diag1 = [midPoints_ver, zeros(obj.m,1)]';
				% Reshape the matrix into a column vector
				diag1 = reshape(diag1, obj.m^2,1);
				% Remove the last element of the column vector to make it a
				% (m^2-1)x1 vector
				diag1 = diag1(1:obj.m^2-1);

				diagM = k("cpx",obj.m);

				% Prepare the vectors to be inserted as diagonals in the sparse
				% matrix.
				d1 = [diagM;zeros(obj.m,1)];
				d2 = [diag1;0];
				d3 = [0;diag1];
				d4 = [zeros(obj.m,1);diagM];
					
				% b is constructed so that every m elements there is the element
				% 2*k_{i,1}
				tempTable = zeros(obj.m,obj.m);
				tempTable(1,:) = k("bl",obj.m);
				b = 2*reshape(tempTable,obj.m^2,1);

				% The main diagonal of A is such that every m elements there is the
				% element 2*k_{i,1}, every element (traslated by 1) there is the
				% element 2*k_{i,m} and then there is the sum of all the other
				% diagonals. 
				tempTable(m,:) =  k("br",obj.m); 
				diag0 = 2*reshape(tempTable,obj.m^2,1)+d1+d2+d3+d4;

				A = spdiags([-d1,-d2,diag0,-d3,-d4],[-obj.m,-1,0,1,obj.m], ...
					obj.m^2,obj.m^2);

				obj.solutionPoints = A\b;
			end
		end

		function value = getSolutionValue(obj,x, varargin)
			if obj.d==1
				if not(all(x>=0)) || not(all(x<=1))
					error("Invalid input: x must be in the range [0,1].");
				end 
				value = obj.solutionPoints(max(1,ceil(x*obj.m)));
			else
				if nargin ~= 3
					error("The function needs to receive in input also the y " + ...
						"values");
				end
				y = varargin{1};
				if any(x<0 | x>1,"all") 
					error("Invalid input: x must be in the range [0,1]");
				end
				if any(y<0 | y>1,"all") 
					error("Invalid input: y must be in the range [0,1]");
				end
				% The first part of epression finds the row, while the second part
				% of the expression finds the column. Given that solutionPoints is
				% vector, the row must be multiplied by m and then we have to sum
				% the column number.
				value = obj.solutionPoints( ...
					(max( 1, ceil(y*obj.m))-1)*obj.m + max(1,ceil(x*obj.m)) );
			end
		end
	end
end

