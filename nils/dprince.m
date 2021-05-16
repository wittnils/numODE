% b_1 and b_2 are the different quadrature weights. 
function y = rungeKutta(f, I, initialvalue, T, epsilon)
    A = [0,0,0,0,0,0,0;
         0.2,0,0,0,0,0,0;
         3./40,9./40,0,0,0,0,0;
         44./45,-56./15,32./9,0,0,0,0;
         19372./6561,-25360./2187,64448./6561,-212./729,0,0,0;
         9017./3168,-355./33,46723./5247,49./176,-5103./18656,0,0;
         35./384,0,500./1113,125./192,-2187./6784,11./84,0]
    b_tilde = [35./384,0,500./1113,125./192,-2187./6784,11./84,0];
    b = [5179./57600,0,7571./16695,393./640,-92097./339200,187./2100,1./40];
    c = [0,1./5,3./10,4./5,8./9,1,1]

    % dimension of the IVP
    d = size(initialvalue)(1)

	% some convienince transpositions
	if size(initialvalue,2) > 1
		initialvalue = initialvalue';
	end

	if size(b,1) > 1
		b = b';
	end

	if size(c,1) > 1
		c = c';
	end

	
    x = zeros(d,2);
	
	y = zeros(d,n);
    y_tilde = zeros(d,n);
	
	k = zeros(d,s);

end

% turn this into a funciton: s = 7 since DP45 is a 7-stage-method
% h=current step size, A = matrix for coeffs, b,b_tilde = quadrature weights, c=interpolation points, (t_0,y_0) initial data
function x = one_step(f,h,A,b,b_tilde,c,d,y0,t0) 
    k= zeros(d,7);
	for i = 1:7
		% compute g_i
		g = y0;
		if i>=2
			% add computed k, scaled with coeff of A
			g += h*sum(k(:,1:i-1).*A(i,1:i-1),2);
		end
		% evaluate f to compute next k_i
		k(:,i) = f(g, t0+c(i)*h);
	end
	% found y(t_j) and y_tilde(t_j)
	y = y0 + h * sum((k.*b),2);
    y_tilde = y0 + h * sum((k.*b_tilde),2);
    x = [y,y_tilde]; 
end
