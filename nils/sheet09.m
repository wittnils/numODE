% Backtracking line search for finding a root of the gradient of the Rosenbrock function F 
% The gradient is given by grad(F)(x,y) = (400x^3-400xy+2x-2,200(y-x^2))^t  
function sheet09()  
    [y1,d1,alpha1] = newton(@f,@jacobian,[1.2;1.2],0.001);
    [y2,d2,alpha2] = newton(@f,@jacobian,[0.001;0.005],0.001); 
    figure 
    hold on 
    plot(y1(1,:),y1(2,:))
    plot(y2(1,:),y2(2,:))
end

function [x,d,alpha] = newton(f,jacobian,x0,epsilon)
    x = [x0];
    d = []; 
    alpha = []; 
    k = 0; 
    while !(norm(f(x(:,end))) <= epsilon)
        k = 0; 
        d(:,end+1) = -inv(jacobian(x(:,end)))*f(x(:,end));
        while (norm(f(x(:,end)+2.^(-k)*d(:,end))) >= norm(f(x(:,end))))
            k = k+1; 
        end 
        alpha(:,end+1) = 2.^(-k); 
        x(:,end+1) = x(:,end) + 2.^(-k)*d(:,end);  
    end
end 

% The jacobian is given by jac(f)(x,y) = (1200x^2-400y+2,-400x \\ -400x,200)
function jac = jacobian(x)
        jac(1,1) = 1200*x(1).^2-400*x(2)+2; 
        jac(1,2) = -400*x(1); 
        jac(2,1) = -400*x(1);
        jac(2,2) = 200; 
end 

% f = grad(F)
function y = f(x) 
    y(1) = 400*x(1).^3-400*x(1)*x(2)+2*x(1)-2; 
    y(2) = 200*(x(2)-x(1).^2);
    y = y'; 
end 