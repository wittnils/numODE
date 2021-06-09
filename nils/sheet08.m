# Given the map H : R^2 \to R, H(p,q) = 1/2p^2+1/2q^2
# We have the system of ODEs
# (p',q') = (-q,p) 
function sheet08(t_0,t_1,x0,y0,h)
    n = abs(t_1-t_0)/h;
    T = linspace(t_0,t_1,n);
    y = implicitEuler(t_0,t_1,x0,y0,h);
    z = explicitEuler(t_0,t_1,x0,y0,h);
    v = symplecticEuler(t_0,t_1,x0,y0,h)
    figure  
    hold on 
    plot3(y(1,:),y(2,:),T,"linewidth",2)
    plot3(z(1,:),z(2,:),T,"linewidth",2)
    plot3(v(1,:),v(2,:),T,"linewidth",2)
    legend("Implicit Euler","Explicit Euler","Symplectic Euler")
    view(3)
    hold off
end

function y = implicitEuler(t_0,T,x0,y0,h)
    n = abs(T-t_0)/h;
    y = [[x0;y0]];
    x_iter = x0;
    y_iter = y0;  
    for i = 1:n-1
        z = [1./(1+h^2)*(x_iter-h*y_iter);y_iter+h*(1./(1+h^2)*x_iter-h*y_iter)]; 
        y(:,end+1)=[z(1),z(2)]; 
        x_iter = z(1);
        y_iter = z(2); 
    end 
end 

function y = explicitEuler(t_0,T,x0,y0,h)
    n = abs(T-t_0)/h;
    y = [[x0;y0]];
    x_iter = x0;
    y_iter = y0;  
    for i = 1:n-1
        z = [x_iter-h*y_iter;y_iter+h*x_iter]; 
        y(:,end+1)=[z(1),z(2)]; 
        x_iter = z(1);
        y_iter = z(2); 
    end 
end

function y = symplecticEuler(t_0,T,x0,y0,h)
    n = abs(T-t_0)/h;
    y = [[x0;y0]];
    x_iter = x0;
    y_iter = y0;  
    for i = 1:n-1
        z = [x_iter-h*y_iter;y_iter+h*(x_iter-h*y_iter)]; 
        y(:,end+1)=[z(1),z(2)]; 
        x_iter = z(1);
        y_iter = z(2); 
    end 
end