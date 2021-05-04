%
% witt_borowitz_sheet2.m
%  
% 
%

function witt_borowitz_sheet2

disp("   k         h     yhN_k    Diff     a_k");

K = 10;
y = zeros(K+1,5);

for k = 0:K % choose upper bound

Nk = 2^(k+1);
yhN = explicitEuler(@(x) cos(pi/2 * x) - 2*x, 0, 1, 1/(2^k), Nk); 

if k == 0
	y(1,:) = [0, 1, yhN, 0, 0];
elseif k==1 
	y(2,:) = [1, 0.5, yhN, yhN-y(1,3), 0];
else
	y(k+1,:) = [k, 1/(2^k), yhN, yhN-y(k,3), (1/log(2)) * log(abs((y(k-1,3)-y(k,3))/(y(k,3)-yhN)))];
end

end

disp(y)

end


function yhN = explicitEuler(f, t0, u0, h, N)
%
% solve autonomous IVP with
% with the explicit Euler method
% 
% @param f function-handle
% @param t0 position of initial value
% @param u0 initial value
% @param h interval size
% @param N steps
%
% @return yhN estimating u(t_0+N*h)

yhN = u0;
for i = 1:N
	yhN = yhN + h*f(yhN);
end

end
