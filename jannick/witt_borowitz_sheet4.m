% sheet 4 - Nils Witt, Jannick Borowitz
function witt_borowitz_sheet4(initialvalue, T, A, b, c)
	% solving lottka volterra
	
	figure
	hold on
	y = rungeKutta(@lotka, initialvalue, T, A, b, c);
	y(:,1:5)
	plot(y(1,:), y(2,:));

	%x = lsode("lotka", [0.5, 1], T);
	%plot(x(:,1), x(:,2));

	xlabel("y_1");
	ylabel("y_2");
end

% lottka volterra
function xdot = lotka(x,t)
a=1;
b=1;
c=1;
d=1;
xdot = zeros(2,1);
xdot(1) = a*x(1)-b*x(1)*x(2);
xdot(2) = c*x(1)*x(2)-d*x(2);
end



