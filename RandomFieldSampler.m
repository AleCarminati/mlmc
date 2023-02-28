classdef RandomFieldSampler
	properties
		d % Number of dimensions
		m_kl % Number of terms of the KL expansion
		sigma2 % variance of exponential covariance function
		lambda % correlation length of exponential covariance function
		w % vectors of the solutions of the transcendental equation
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
		function obj = RandomFieldSampler(m_kl, sigma2, lambda)
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
			obj.A_n = sqrt(1./( 0.5 - sin(2.*obj.w)./(4.*obj.w) + ...
				obj.lambda^2/2.*obj.w.*(obj.w+sin(obj.w).*cos(obj.w))+ ...
				obj.lambda.*sin(obj.w).^2 ));

			if d==2
				[obj.idx_1_2d,obj.idx_2_2d]=meshgrid(1:obj.m_kl);
				obj.idx_1_2d = reshape(obj.idx_1_2d,[],1);
				obj.idx_2_2d = reshape(obj.idx_2_2d,[],1);
				prodMatr = tril(obj.eigenvalues*obj.eigenvalues');
				prodVec = reshape(prodMatr,[],1);
				[obj.eigenvalues, idx] = sort(prodVec,'descend');
				% Remove all the zeros of the triangular matrix
				obj.eigenvalues = reshape(obj.eigenvalues,obj.m_kl,1);
				obj.idx_1_2d = reshape(obj.idx_1_2d(idx),m_kl,1);
				obj.idx_2_2d = reshape(obj.idx_2_2d(idx),m_kl,1);
			end
			
			obj = obj.updateRandom();

			obj.b_n = @(x) obj.A_n.*(sin(obj.w.*x)+ ...
				obj.lambda.*obj.w.*cos(obj.w.*x));
		end
		
		function w = findSolutionTransEq(obj)
			% This function finds the first m_kl solutions of a transcendental 
			% equation of the form tan(x)-(2lambda)/(lambda^2x^2+1)=0, using the 
			% Matlab function fzero. The solutions are stored in the vector w.
			% Output:
			% - w: a vector of size m_kl containing the computed solutions

			f = @(x) tan(x)-(2*obj.lambda)./(obj.lambda^2.*x.^2+1);
			w = zeros(obj.m_kl,1);
			idx = 1;
			w(idx) = 0;
			idx=idx+1;
			while idx<=obj.m_kl
				zero_1 = fzero(f,[-pi/2*(idx+1), -pi/2*(idx-1)]);
				zero_2 = fzero(f,[pi/2*(idx-1), pi/2*(idx+1)]);
				w(idx) = min(zero_1,zero_2);
				idx=idx+1;
				if(idx<=obj.m_kl)
					w(idx)=max(zero_1,zero_2);
					idx = idx+1;
				end
			end
		end

		function value = computeRandomFieldValue(obj,x,varargin)
			if obj.d==2 && nargin~=3
				error("Input error. When working in 2D, the input requires " + ...
					"also the second coordinate.");
			end

			% Check that x is a column vector.
			if size(x,2)~=1
				error("Input error. The input of the function to compute " + ...
					"the random field value must be a column vector.");
			end

			y = varargin{1};
			if obj.d==2 
				% Check that y is a column vector.
				if size(y,2)~=1
				error("Input error. The input of the function to compute " + ...
					"the random field value must be a column vector.");
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
		end

		function obj = updateRandom(obj)
			% This function updates the random part of the KL expansion,
			% generating a new random field.
			obj.xi =randn(obj.m_kl,1);
		end
	end
end
