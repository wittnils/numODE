% sheet 4 - Nils Witt, Jannick Borowitz

% README
% initialvalue must be a column vector
% T must be a row vector
% A is a lower triangular matrix
% b,c must be row vectors
% y ist a matrix; every column vector of y represents a discrete solution
%
% EXAMPLE - ERK of order 4
%
%	A = [0,0,0,0;0.5,0,0,0;0,0.5,0,0;0,0,1,0]
%	c = [0,0.5,0.5,1]
%	b = [1/6,2/6,2/6,1/6]
%
%	initialvalue = [0.5;1]
%	T = linspace(0,300,2000);
%
%	witt_borowitz_sheet4(initialvalue, T, A, b, c)
%

function witt_borowitz_sheet4(initialvalue,h)
	% solving lottka volterra
    T = linspace(0,100,100/h);
    A = [0,0,0,0;0.5,0,0,0;0,0.5,0,0;0,0,1,0]
	c = [0,0.5,0.5,1];
    b = [1/6,2/6,2/6,1/6];

	figure
	hold on
	y = rungeKutta(@lorenz, initialvalue, T, A, b, c);
	plot(y(1,:), y(3,:));

	% TESTING against lsode (octave)
	% x = lsode("lotka", initialvalue, T);
	% plot(x(:,1), x(:,2));

	xlabel("x(t)");
	ylabel("z(t)");
	title("Lorentz");
end

% lottka volterra
% implementation from sheet 1
function xdot = lorenz(x,t)
sigma=4;
rho=-1;
bet=0.5;
xdot = zeros(3,1);
xdot(1) = sigma*(x(2)-x(3));
xdot(2) = x(1)*(rho-x(3))-x(2);
xdot(3) = x(1)*x(2)-bet*x(3); 
end

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
% @return y		dxn-sized matrix containing the solution y(:,i)\cong u(t_i)
function y = rungeKutta(f, initialvalue, T, A, b, c)

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

	n = size(T)(2)
	d = size(initialvalue)(1)
	
	y = zeros(d,n);
	
	s = size(c)(2)
	k = zeros(d,s);

	t0 = T(1);
	y(:,1) = y0 = initialvalue;

	% determine every y(t_j)
	for j = 2:n
		t = T(j);
		h = t - t0;
		
		% determine every k_i for y_j
		k = zeros(d,s);
		for i = 1:s
			% compute g_i
			g = y0;
			if i>=2
				% add computed k, scaled with coeff of A
				g += h*sum(k(:,1:i-1).*A(i,1:i-1),2);
			end
			% evaluate f to compute next k_i
			k(:,i) = f(g, t0+c(i)*h);
		end
		% found y(t_j)
		y(:,j) = y0 + h * sum((k.*b),2);
		
		% used for next approx
		t0 = t;
		y0 = y(:,j);
	end

end
