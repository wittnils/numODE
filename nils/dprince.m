% b_1 and b_2 are the different quadrature weights. 
function dprince(I,initialvalue,epsilon)
	% solving lottka volterra
	
	figure
	hold on
	[T,y] = stepsizecontrolled(I,initialvalue,epsilon);
	plot(y(1,:), y(2,:));

	% TESTING against lsode (octave)
	% x = lsode("lotka", initialvalue, T);
	% plot(x(:,1), x(:,2));

	xlabel("y_1");
	ylabel("y_2");
	title("Lotka-Volterra");
end

function [T,y] = stepsizecontrolled(I, initialvalue, epsilon)
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

    T=[I(1)];
    d = size(initialvalue)(1)
	% storing solutions
    x = zeros(d,2);	
	y = [initialvalue];
    % setting initial parameters 
    t = T(1);
    h_optimal = 0;
    % h_opt is = 0, this way we can ensure the algorithm does indeed start. 
    h = min(0.1, I(2)-I(1));;
    j=2; 
    while t < I(2)
        % compute two soln of the IVP w/ given starting values (t,y_0), x(1) has been computed w/ b_tilde and x(2) w/ b
        
        [x1,x2] = one_step(@f,h,A,b,b_tilde,c,d,y(:,j-1),t);
        h_optimal =  h_opt(h,x1,x2,epsilon);

        if norm(x1-x2) > epsilon
            % reject
            h = h_optimal;
        else
            % accept
            if I(2)-t <= 1.1*h_optimal
                h = min(I(2)-t,h_optimal); 
            else 
                h = h_optimal;
            end
            T(j)=t+h;
            y(:,j)=x2;

            t = T(j);
            j+=1;
        end
        % setting h = h_opt since at this point the value for h has not been rejected. 
    end
end

% turn this into a funciton: s = 7 since DP45 is a 7-stage-method
% h=current step size, A = matrix for coeffs, b,b_tilde = quadrature weights, c=interpolation points, (t_0,y_0) initial data
function [x1,x2] = one_step(f,h,A,b,b_tilde,c,d,y0,t0) 
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
	x1 = y0 + h * sum((k.*b),2);
    x2 = y0 + h * sum((k.*b_tilde),2); 
end

function xdot = f(x,t)
a=1;
b=1;
c=1;
d=1;
xdot = zeros(2,1);
xdot(1) = a*x(1)-b*x(1)*x(2);
xdot(2) = c*x(1)*x(2)-d*x(2);
end

function h_optimal = h_opt(h,x1,x2,epsilon)
    h_optimal = h*(epsilon/norm(x1-x2))^(1./8); 
end 