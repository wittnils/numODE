function ExplicitEulerIter(n)
    a = zeros(n,5);
    if n >= 2
        for k = 0:n
            a(k+1,1) = k; 
            a(k+1,2) = 2^(-k); 
            a(k+1,3) = explicitEuler(0,1,a(k+1,2),2^(k+1));
            if k+1 > 1
                a(k+1,4) = a(k+1,3)-a(k,3);
            end
            if k+1 > 2 
                a(k+1,5) = (1./log(2))*log(abs(a(k,4)./a(k+1,4)));
            end
        end
    end
    disp(["In dieser Reihenfolge werden die Werte ausgegeben:" newline
          "k   h   y_k   y_k-y_{k-1}   a_k" newline])
    disp(a)
end

function y_n = explicitEuler(t_0,u_0,h,N)
    y_n = u_0; 
    for n = 1:N 
        y_n = y_n + h*(cos(y_n*pi/2)-2*y_n);
    end
end  
