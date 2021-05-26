function sheet06(t0,tmax,u0,h)

	A = [-21 19 -20; 19 -21 20; 40 -40 -40]
    
    	[T,y] = implicitEuler(A,t0,tmax,u0,h);
    	[T2,x] = explicitEuler(A,t0,tmax,u0,h);
    
    	N = (tmax-t0)/h;
	U = zeros(3,N+1);
    	for t = 0:N
		U(:,t+1) = u(t*h+t0);
	end
	disp("time steps")
	[t0:h:tmax]
	disp("exact numerical solutions U")
	U
	disp("implicit euler solutions y")
    	y
	disp("explicit euler solutions x")
    	x
	disp("error abs(y-U)")
	abs(y-U)
	disp("error abs(x-U)")
	abs(x-U)

    #figure
    #hold on
    #t = [t0:h:tmax]
    #plot(t, U)

end

% t_0 is inital time value, T is end time point, u_0 is inital value as column vector, h is step size, A is square matrix of the linear IVP 
function [T,y] = implicitEuler(A,t_0,T,u0,h)
    n = abs(T-t_0)/h;
    T = linspace(t_0,T,n);
    d = size(u0,1);  
    I = eye(d);
    M = inv(I-h*A); 
    y = [u0];
    y_n = y(:,1);  
    for i = 1:n
        y_n = M*y_n; 
        y(:,end+1)=y_n; 
    end 
end 

% Modified & Explicit Euler
function [T,y] = explicitEuler(A,t_0,T,u0,h)
    n = abs(T-t_0)/h;
    T = linspace(t_0,T,n);
    d = size(u0,1); 
    y = [u0];
    y_n = y(:,1);
    for n = 1:n
        y_n = y_n + h*A*(y_n+h/2*A*(y_n));
        y(:,end+1) = y_n;
    end 
end 

function U = u(t)
	U = [0.5 * exp(-2*t) + 0.5 * exp(-40*t) * (cos(40*t) + sin(40*t));0.5 * exp(-2*t) - 0.5 * exp(-40*t) * (cos(40*t) + sin(40*t));-exp(-40*t) * (cos(40*t) - sin(40*t))];	
end
