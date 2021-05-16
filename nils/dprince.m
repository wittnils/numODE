% b_1 and b_2 are the different quadrature weights. 
function [T,y] = rungeKutta(f, I, initialvalue, epsilon)
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
    
    
    % dimension of the IVP, starting value t_0 is left end of intervall I
    d = size(initialvalue)(1)
    T=[I(1)];

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

	% storing solutions
    x = zeros(d,2);	
	y = [initialvalue];

    % setting initial parameters
    y_0 = initialvalue; 
    t = T(1);

    % h_opt is = 0, this way we can ensure the algorithm does indeed start. 
    h = 0.1;

    while t < I(1) 
        % compute two soln of the IVP w/ given starting values (t,y_0), x(1) has been computed w/ b_tilde and x(2) w/ b
        x = one_step(f,h,A,b,b_tilde,c,d,y_0,t);
        while abs(x(1)-x(2)) > epsilon
            h = h*(epsilon/abs(x(1)-x⁽(2)))^(1./8);
            x = one_step(f,h,A,b,b_tilde,c,d,y_0,t);
        end
        
        % setting h = h_opt since at this point the value for h has not been rejected. 
        if I(1)-t <= 1.1*h*(epsilon/abs(x(1)-x⁽(2)))^(1./8)
            h = I(1)-t; 
        end
        h = h*(epsilon/abs(x(1)-x⁽(2)))^(1./8);

        % store time value
        T = [T;t+h];
        % add y_tilde to the soln vector 
        y = [y;x(1)];

        % set parameters for next computation
        t = T(rows(T));
        y0 = y(rows(y));
    end
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
    x = [y_tilde,y]; 
end

function y = f(x,t)
    y = cos(x/2)-t;
end 