% Sheet 5 - Witt, Borowitz
function [T, y] = witt_borowitz_sheet5(f, I, u0, eps)


	%output_precision(10)

	A = 	[0,0,0,0,0,0,0;
		1/5,0,0,0,0,0,0;
		3/40,9/40, 0,0,0,0,0;
		44/45, -56/15, 32/9, 0, 0, 0, 0;
		19372/6561, -25360/2187, 64448/6561, -212/729, 0,0,0;
		9017/3168, -355/33, 46732/5247, 49/176, -5103/18656, 0 ,0;
		35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0]
	b1 = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0]
	b2 = [5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40]
	c = [0, 1/5, 3/10, 4/5, 8/9, 1, 1]
	
	
	eps = 0.000001	

	[T,y] = addaptiveStepSizeControlApprox(@f, I, u0, eps, A, b1, b2, c, 4);
	
	if size(u0, 2) > 1
		y = y';
	end
end

function [T,y] = addaptiveStepSizeControlApprox(f, I, u0, eps, A, b1, b2, c, p)


	y = [];
	T = [];

	if size(u0,2) > 1
		u0 = u0';
	end

	if size(b1,1) > 1
		b1 = b1';
		b2 = b2';
	end

	if size(c,1) > 1
		c = c';
	end

	% b1 has higher consistency order than b2 (p)

	% det h_op, y1 is the sol of the better consistency order
	HOpt = @(y1, y2, h) h*((eps / norm(y2-y1))^(1/(p+1)));

	h = min(0.1, I(2)-I(1));
	t = I(1);
	i = 2;
	y(:,1) = u0;
	T(1) = t;

	while t < I(2) && t+h!=t
		%1
		[y1, y2] = embeddedRungeKutta(f, y(:,i-1), t, h, A, b1, b2, c);
		%2
		hOpt = HOpt(y1,y2,h);
		if isinf(hOpt)
			disp("hOpt is Inf")
			%break;
		end
		%3
		if norm(y2-y1) > eps
			% reject
			h = hOpt;
		else
			% accept
			y(:,i) = y1;
			T(i) = t+h;
			t = T(i);

			if I(2)-t <= 1.1 * hOpt || isinf(hOpt)
				h = I(2)-t;
			elseif y1 != y2
				h = hOpt;
			end
			% 5
			i += 1;
		end
	end

end

function [y1, y2] = embeddedRungeKutta(f, initialvalue, t, h, A, b1, b2, c)
	d = size(initialvalue)(1);

	s = size(c)(2);
	k = zeros(d,s);

	y0 = initialvalue;
		
	k = zeros(d,s);
	for i = 1:s
		g = y0;
		if i>=2
			g += h*sum(k(:,1:i-1).*A(i,1:i-1),2);
		end
		k(:,i) = f(g, t+c(i)*h);
	end
	y1 = y0 + h * sum((k.*b1),2);
	y2 = y0 + h * sum((k.*b2),2);
		
end
