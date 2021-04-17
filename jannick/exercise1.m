%  Exercise 1
%  Solving Lotka-Volterra equations:
%	u' 	= a*u - b*u*v
%	v' 	= b*u*v - c*v 
%  with a=b=c=d=1 and starting values
%  	u(0) 	= 0.5*j,	j=1,...,10
%	v(0) 	= 1
%
function exercise1
% discrete equidistant time steps
t = linspace(0, 300, 1000);

% start plotting
figure
hold on

% determine approx
for i = 1:10
	x = lsode("lotka", [0.5*i, 1], t);
	plot(x(:,1), x(:,2));
end

% plot configurations
xlabel("x_1");
ylabel("x_2");
end

function xdot= lotka(x,t)
a=1;
b=1;
c=1;
d=1;
xdot = zeros(2,1);
xdot(1) = a*x(1)-b*x(1)*x(2);
xdot(2) = c*x(1)*x(2)-d*x(2);
end
