

function test2_rungekutta_sheet4

	A = [0,0,0,0;0.5,0,0,0;0,0.5,0,0;0,0,1,0]
	c = [0;0.5;0.5;1]
	b = [1/6;2/6;2/6;1/6]

	T = linspace(0,300,1000);

	for i = [1:3]
		initialvalue = [0.5*i,1]
		witt_borowitz_sheet4(initialvalue, T, A, b, c)
	end
end
