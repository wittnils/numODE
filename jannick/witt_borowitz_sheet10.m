function witt_borowitz_sheet10()
	%output_precision(10);
	f = @(x, t) -200*t*x^2
	u0 = 1/901
	t0 = -3
	% exact solution
	u = @(t) 1./(1+100*t.^2)
	t = 0
	H = 1./2.^[0:8]
	U = u(0)

	disp_results(U, @lmmNaive, f, u0, t0, t, H)
	disp_results(U, @lmmSkeem, f, u0, t0, t, H)
	disp_results(U, @rk, f, u0, t0, t, H)

end

% runge kutta of appropriate order
% approx steps [t0:h:t]
function [T, Y, c] = rk(f, h, t, t0, u0)
	A = [0,0,0,0;0.5,0,0,0;0,0.5,0,0;0,0,1,0];
	c = [0,0.5,0.5,1];
	b = [1/6,2/6,2/6,1/6];	

	T = [t0:h:t];
	Y = rungeKutta(f, [u0], T, A,b,c);
	% 4 evaulations in one step, since rk is an 4 stage rk
	c = 4*(size(T,2)-1);
end

function [T, y, c] = lmmSkeem(f, h, t, t0, u0)
	T = [t0:h:t];
	alpha = [0,0,0,-1];
	beta = [-9, 37, -59, 55]*(h/24);
	s = [0,0,0,0];
	
	% our update handle
	update_s = @(ti, yi, s2) [0, s2(1:3)] + alpha*yi - beta*f(yi, ti);
	
	y = [];
	%s = update_s(t0, u0, s);
	% use RK for 3 approx
	[t1, y1, d] = rk(f, h, T(4), t0, u0);
	c = d;
	for i = [1:4]
		s = update_s(T(i), y1(i), s);
		y(end+1) = y1(i);
	end
	for ti = T(5:end);
		yi = -s(end);
		s = update_s(ti, yi, s);
		y(end+1) = yi;
		c += 1;
	end
end

function [T, y, c] = lmmNaive(f, h, t, t0, u0)
	T = [t0:h:t];
	[t1, y1, d] = rk(f, h, T(4), t0, u0);
	c = d;
	F = [];
	for i = [1:4]
		F(end+1) = f(y1(i), T(i));
	end
	y = y1';
	for ti = T(5:end)
		y(end+1) = y(end) + h/24 * sum([-9, 37, -59, 55].*F);
		F = [F(2:4), f(y(end),ti)];
		c += 1;
	end
	y = y';
end
% disp results
% initial value u(t0) = u0
function disp_results(U, m, f, u0, t0, t, H)

	disp(m)
	disp(list_in_columns({"   -----  h  -----"," -----error-----","-----evals-----"}))
	res = [];
	for i = [1:size(H,2)]
		[T, Y, c] = m(f, H(i), t, t0, u0);
		res(end+1,:) = [H(i), abs(Y(end)-U), c];
	end
	disp(res)
end
