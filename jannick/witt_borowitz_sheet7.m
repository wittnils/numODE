function witt_borowitz_sheet7

	A1 = [1, 2; 1, 0];
	A2 = [1, 2; -2, 1];
	A3 = [-1 20; -20, -1];

	a = 5
	for u0 = 1:a
		for u1 = 1:a
			initial_values(u0*a+u1,:)=[u0-2,u1-2];
		end
	end
	
	plot_trajectories(initial_values, A1, "A1");
	plot_trajectories(initial_values, A2, "A2");
	plot_trajectories(initial_values, A3, "A3");
end

function plot_trajectories(initial_values, A, figure_title)
	n = size(initial_values, 1)
	N = 300
	T = linspace(1,0.1,N);

	figure
	hold on
	for i = 1:n
		[ts, Y] = ode45(@(t, y) A*y, T, initial_values(i,:));
		plot(Y(:,1), Y(:,2));
	end
	title(figure_title)
	hold off
end
