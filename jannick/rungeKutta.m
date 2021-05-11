% solving ODE using abitrary ODE
% @param f		function f(x,t) RHS of ODE: u'(t)=f(u(t),t)
% @param initialvalue	d-sized vector containing initial properties for t0 = 0 
% @param T		n-sized vector containing strictly inceasing time steps t1,..,tn
% @param A		s x s lower triangular matrix (for any s)
% @param b		s-sized vector
% @param c		s-sized vector
%
% A,b,c represent the Butcher tableau
% 
% @return y		n-sized vector containing the solution y(i)\cong u(t_i)
function y = rungeKutta(f, initialvalue, T, A, b, c)
	n = size(T)(2)
	d = size(initialvalue)(1)
	y = zeros(d,n);
	
	s = size(c)(2)
	k = zeros(d,s);

	t0 = T(1);
	y(:,1) = y0 = initialvalue;
	for j = 2:n
		t = T(j);
		h = t - t0;
		
		k = zeros(d,s);
		for i = 1:s
			g = y0;
			if i>=2
				g += h*sum(k(:,1:i-1).*A(i,1:i-1),2);
			end
			k(:,i) = f(g, t0+c(i)*h);
		end
		y(:,j) = y0 + h * sum((k.*b),2);
		
		t0 = t;
		y0 = y(:,j);
	end

	%y(:,1:5)
end
