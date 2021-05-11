%source (rungeKutta)

function test_rungekutta_sheet4

	A = [0,0;0.5,0];
	b = [0,1];
	c = [0,0.5];


	f = @(y, t) cos((pi/2) * y) - 2*y;

	u0 = [1];

	%y = zeros(1,10);
	y = rungeKutta(f, u0, [0:0.01:2], A, b, c);

	disp("res");
	y
end

%
% Modified Euler - One Step Method
%
function yN = modifiedEuler(f, t0, u0, h, N)

	yN = u0;
	k1 = 0;
	k2 = 0;

	while N>0
		k1 = f(t0, yN);
		k2 = f(t0+0.5*h, yN+0.5*h*k1);
		yN = yN + h*k2;
		
		N-=1;
	end
end


%
% printAnlysis - print Analysis for one step method 
% that approximates an sol of an IVP at u(T)
%
% @param osm 	function handle that satifies osm(f, t0, u0, h, N)
% @param t0 	initial position
% @param t1 	initial value
% @param K 	K+1 refinements of h (time-step-size)
% @param T	T position to approximate
%
function printAnalysis(osm, f, t0, u0, K, T)

	% matrix which stores result
	R = zeros(K+1, 5);

	% lambda function for caculating convergence rate
	a = @(k, y, R) (1/log(2)) * log(abs((R(k-1,3)-R(k,3))/(R(k,3)-y)));

	% compute for every h=2^(-k)
	for k = 0:K
		h = 2^(-k);
		N = T/h;
		yN = osm(f, t0, u0, h, N);

		if k==0
			R(1,:) = [k, h, yN, 0, 0];
		elseif k==1
			R(2,:) = [k, h, yN, yN-R(1,3), 0];
		else
			R(k+1,:) = [k, h, yN, yN-R(k,3), a(k, yN, R)];
		end
	end

	disp("   k         h        yhN_k       Diff       a_k");
	disp(R);

end
