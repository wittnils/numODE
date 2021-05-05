function sheet02(n)
    % a ~ Modified Euler, b ~ Heun 
    a = zeros(n,5);
    b = zeros(n,5);
    if n >= 2
        for k = 0:n
            a(k+1,1) = k; 
            a(k+1,2) = 2^(-k); 
            b(k+1,1) = k; 
            b(k+1,2) = 2^(-k); 
            a(k+1,3) = modEuler(0,1,a(k+1,2),@f);
            b(k+1,3) = heun(0,1,b(k+1,2),@f);
            if k+1 > 1
                a(k+1,4) = a(k+1,3)-a(k,3);
                b(k+1,4) = b(k+1,3)-b(k,3);
            end
            if k+1 > 2 
                a(k+1,5) = (1./log(2))*log(abs(a(k,4)./a(k+1,4)));
                b(k+1,5) = (1./log(2))*log(abs(b(k,4)./b(k+1,4)));
            end
        end
    end
    disp(["In dieser Reihenfolge werden die Werte ausgegeben:" newline
          "k   h   y_k   y_k-y_{k-1}   a_k" newline])
    a
    disp(["In dieser Reihenfolge werden die Werte ausgegeben:" newline
          "k   h   y_k   y_k-y_{k-1}   b_k" newline])
    b
end

function y = f(x)
    y = cos(x*pi/2)-2*x;
end

% Heun-Method
function y_n = heun(t_0,u_0,h,f)
    y_n = u_0; 
    N = 2./h;
    k_1 = 0;
    k_2 = 0;
    for n = 1:N 
        % In the second summand we use that f is independent of t, thus we omit it in the function parameter list 
        y_n =y_n+h*1./2*(f(y_n)+f(y_n+h*f(y_n)));
    end
end 

% Modified Euler
function y_n = modEuler(t_0,u_0,h,f)
    y_n = u_0;
    N = 2./h; 
    k_1 = 0;
    k_2 =0 ;
    for n = 1:N
        y_n = y_n + h*f(y_n+f(y_n)*h/2);
    end 
end 