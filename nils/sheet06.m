function [y,x] = sheet06(A,t0,tmax,u0,h)
    figure 
    hold on 
    [T,y] = implicitEuler(A,t0,tmax,u0,h);
    [T2,x] = explicitEuler(A,t0,tmax,u0,h);
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
    for i = 1:n-1
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
    for n = 1:n-1
        y_n = y_n + h*A*(y_n+h/2*A*(y_n));
        y(:,end+1) = y_n;
    end 
end 

