# Given the map H : R^2 \to R, H(p,q) = 1/2p^2+1/2q^2
# We have the system of ODEs
# (p',q') = (-q,p) 
function sheet08(t_0,T,x0,y0,h)
    [T,y] = implicitEuler(t_0,T,x0,y0,h);
    figure  
    hold on 
    plot(y(1,:),y(2,:))
    hold off
end

function [T,y] = implicitEuler(t_0,T,x0,y0,h)
    n = abs(T-t_0)/h;
    T = linspace(t_0,T,n);
    y = [[x0;y0]];
    x_iter = x0;
    y_iter = y0;  
    for i = 1:n
        z = [1./(1+h^2)*(x_iter-h*y_iter);y_iter+h*(1./(1+h^2)*x_iter-h*y_iter)]; 
        y(:,end+1)=[z(1),z(2)]; 
        x_iter = z(1);
        y_iter = z(2); 
    end 
end 