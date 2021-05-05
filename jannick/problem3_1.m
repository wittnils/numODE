% problem3_1
function problem3_1

	R = zeros(7,6);

	a = @(h,i) h*s(i);
	convg1 = @(h, i, a) (1/log(2))*log(abs(a(h,i)/a(h/2,i+1)));
	
	convg2 = @(h, i, a) (1/log(2))*log(abs((a(h,i) - a(h/2,i+1))/(a(h/2,i+1)-a(h/4,i+2))));

	for i = 0:6
		h = 1/(2^i);
		R(i+1,:) = [i, h, a(h,i), 0, convg1(h, i, a), convg2(h, i, a)];
		if i > 0
			%R(i+1,4)
			v = R(i+1,:);
			v(4) = abs(a(h,i)-a(2*h,i-1));
			R(i+1,:) = v;		
		end
	end
	disp(R);

end

function y = s(i)

	if mod(i,4) == 0 || mod(i,4) == 3
		y = 1;
	else
		y = -1;
	end

	y
end
