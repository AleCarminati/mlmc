classdef RandomFieldSampler
	% It implements a random field generator using the Karhunen-Loeve 
	% expansion (KL expansion) with exponential covariance function.

	properties
		d % Number of dimensions (1 or 2)
		m_kl % Number of terms of the KL expansion
		sigma2 % variance of the exponential covariance function
		lambda % correlation length of the exponential covariance function
		w % vector of the solutions of the transcendental equation
		eigenvalues % vector of the eigenvectors of the covariance function
		xi % vector of the independent standard Gaussian variables of the KL 
			% expansion 
		A_n % vector of the normalizing coefficents of the eigenfunctions of 
			% the KL expansion in 1D
		b_n % vector of the eigenfunctions of the KL expansion in 1D
		idx_1_2d % Indexes of the first 1D eigenvalue that is used in the 
			% product to obtain the 2D eigenvalue
		idx_2_2d % Indexes of the second 1D eigenvalue that is used in the 
			% product to obtain the 2D eigenvalue
	end
	
	methods
		function obj = RandomFieldSampler(m_kl, sigma2, lambda, d)
			if d~=1 && d~=2
				error("Invalid input: d=%d, but d must be equal to 1 or 2.", d);
			end
			if sigma2<0
				error("Invalid input: sigma^2= %d, but the variance of the " + ...
					"exponential covariance function must be positive.");
			end
			obj.d = d;
			obj.m_kl = m_kl;
			obj.sigma2 = sigma2;
			obj.lambda = lambda;

			obj.w = obj.findSolutionTransEq();
			obj.eigenvalues = 2*obj.lambda./(obj.lambda^2.*obj.w.^2+1);
			obj.A_n = sqrt(4.*obj.w./ (...
				2.*obj.w.*(obj.lambda+obj.lambda.^2.*obj.w.^2+1) + ...
				(obj.lambda.^2.*obj.w.^2-1).*sin(2.*obj.w)- ...
				2.*obj.lambda.*obj.w.*cos(2.*obj.w) ...
			));

			if d==2
				[obj.idx_1_2d,obj.idx_2_2d]=meshgrid(1:obj.m_kl);
				obj.idx_1_2d = reshape(obj.idx_1_2d,[],1);
				obj.idx_2_2d = reshape(obj.idx_2_2d,[],1);
				prodMatr = obj.eigenvalues*obj.eigenvalues';
				prodVec = reshape(prodMatr,[],1);
				[obj.eigenvalues, idx] = sort(prodVec,'descend');
				% Remove all the zeros of the triangular matrix
				obj.eigenvalues = obj.eigenvalues(1:obj.m_kl);
				obj.idx_1_2d = obj.idx_1_2d(idx);
				obj.idx_1_2d = obj.idx_1_2d(1:obj.m_kl);
				obj.idx_2_2d = obj.idx_2_2d(idx);
				obj.idx_2_2d = obj.idx_2_2d(1:obj.m_kl);
			end

			obj.eigenvalues = obj.sigma2.*obj.eigenvalues;
			
			obj = obj.updateRandom();

			obj.b_n = @(x) obj.A_n.*(sin(obj.w.*x)+ ...
				obj.lambda.*obj.w.*cos(obj.w.*x));
		end
		
		function value = computeRandomFieldValue(obj,x,varargin)
			if obj.d==2 && nargin~=3
				error("Input error. When working in 2D, the input requires " + ...
					"also the second coordinate.");
			end

			% Check that x is a column vector or a scalar.
			if not(isscalar(x)) && size(x,2)~=1
				error("Input error. The input of the function to compute " + ...
					"the random field value must be a column vector or a scalar.");
			end

			if obj.d==2 
				y = varargin{1};
				% Check that y is a column vector or a scalar.
				if not(isscalar(y)) && size(y,2)~=1
				error("Input error. The input of the function to compute " + ...
					"the random field value must be a column vector or a scalar.");
				end
				% Check that x and y are compatible.
				if length(x)~=length(y)
					error("Input error. The inputs of the function to compute" + ...
						" the random field value must have the same length");
				end
			end

			value = zeros(length(x),1);
			if obj.d == 1
				for i=1:length(x)
					value(i) = sum(obj.xi.*sqrt(obj.eigenvalues).*obj.b_n(x(i)));
				end
			else
				for i=1:length(x)
					eigfun_res_x = obj.b_n(x(i));
					eigfun_res_x = eigfun_res_x(obj.idx_1_2d);
					eigfun_res_y = obj.b_n(y(i));
					eigfun_res_y = eigfun_res_y(obj.idx_2_2d);
					value(i) = sum(obj.xi.*sqrt(obj.eigenvalues(obj.idx_1_2d) ...
						.*obj.eigenvalues(obj.idx_2_2d)).*eigfun_res_x.*eigfun_res_y);
				end
			end
			value = exp(value);
		end

		function obj = updateRandom(obj)
			% This function updates the random part of the KL expansion,
			% generating a new random field.
			obj.xi = randn(obj.m_kl,1);
		end
	end

	methods(Access=private)
		function w = findSolutionTransEq(obj)
			% This function finds the solutions to a transcendental equation 
			% using the fzero function in Matlab.
			% Output:
			% - w: a vector of size m_kl containing the computed solutions

			f = @(x) tan(x)-(2*obj.lambda.*x)./(obj.lambda^2.*x.^2-1);
			w = zeros(obj.m_kl,1);
			% To avoid problems with tan(x) near pi/2(idx+1), we add a little
			% epsilon.
			eps = 1e-12;
			idx_increased = false;
			n_interval = 1;
			if 1/obj.lambda<pi/2
				% The exceptional case when the interval (0,pi/2] contains a zero. 
				w(n_interval) = fzero(f,[1/obj.lambda+eps, pi/2-eps]);
				idx_increased = true;
			end
			while n_interval+idx_increased<=obj.m_kl
				if pi/2*(2*n_interval-1)+eps<1/obj.lambda && ...
					1/obj.lambda<pi/2*(2*n_interval+1)-eps
					% The exceptional case when the interval contains two zeros.
					w(n_interval) = fzero( ...
						f,[pi/2*(2*n_interval-1)+eps, 1/obj.lambda-eps]);
					w(n_interval+1) = fzero( ...
						f,[1/obj.lambda+eps, pi/2*(2*n_interval+1)-eps]);
					n_interval = n_interval + 1;
					idx_increased = true;
				else
					% Since f is an odd function and the eigenvalues are computed as 
					% the square of the solution, only the positive solutions
					% are returned.  
					w(n_interval+idx_increased) = fzero( ...
						f,[pi/2*(2*n_interval-1)+eps, pi/2*(2*n_interval+1)-eps]);
					n_interval = n_interval+1;
				end
			end
		end
	end
end
