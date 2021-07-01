function witt_borowitz_sheet11

	coeff1 = @(z) [1, -1-z];
	coeff2 = @(z) [1, -1-1.5*z, 0.5*z];
	coeff3 = @(z) [1, -1-(23/12)*z, (4/3)*z, (-5/12)*z];

	plotStabilityRegion(coeff1)
	plotStabilityRegion(coeff2)
	plotStabilityRegion(coeff3)
end


function plotStabilityRegion(coeff)
	z = [];
	n = 3;
	density = 0.05;
	figure
	for x = [-n:density:n]
		for y = [-n:density:n]
			a = complex(x,y);
			if rootTest(roots(coeff(a)))==1
				z(end+1) = a;
			end
		end
	end
	#plot(real(z), imag(z), 'linewidth', 2)
	scatter(real(z), imag(z));
	title("Stability Region")
	axis([-n n -n n]);
	grid on;
end


function ans = rootTest(nullstellen)

	ans = 1;
	for i = 1:size(nullstellen,1)
		if abs(nullstellen(i)) > 1
			ans = 0;
			return
		end
		if abs(nullstellen(i)) == 1
			for j = i+1:size(nullstellen,1)
				if nullstellen(j) == nullstellen(i)
					ans=0;
					return
				end
			end
		end
	end
end
