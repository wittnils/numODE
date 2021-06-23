function sheet10()
    h = 1./2^4;
    T = linspace(-3,0,3./h);
    y1 = rungeKutta(@f,1./901,T);
    y2 = AdamsBashford(@f,T,1./901,h);
    y3 = lsode("f",1./901,T);
    hold on 
    plot(T,y2,"linewidth",2);
    plot(T,y1,"linewidth",2);
    plot(T,y3,"linewidth",2); 
    legend("RKM","ABM","lsode")
    hold off 
    y2
end 

% (t0,u0) initial value, h step size
function y = AdamsBashford(f,T,u0,h)
    % computing first four steps using RKM
    T2 = linspace(-3,-3+4*h,4);
    y = rungeKutta(f,u0,T2); 
    % using ABM to compute the remaining values
    i = 4; 
    t = T(1)+4*h; 
    for i = 1:length(T)-4
        y(end+1) = y(end)+h/24*(55*f(y(end),t)-59*f(y(end-1),t-h)+37*f(y(end-2),t-2*h)-9*f(y(end-3),t-3*h)); 
        t = t+h; 
    end
end

function y = rungeKutta(f, initialvalue, T)
    % butcher-tableau of classical RKM
    A = [0,0,0,0;
        0.5,0,0,0;
        0,0.5,0,0;
        0,0,1,0];
    b = [1./6; 1./3;1./3;1./6]; 
    c = [0;0.5;0.5;1];

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

function y = f(x,t)
    y = -200*t*x^2;
end