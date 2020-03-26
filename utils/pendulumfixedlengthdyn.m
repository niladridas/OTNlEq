function dX = pendulumfixedlengthdyn(~,X,params)
	g = params.g;
	L = params.L;
    dX = zeros(4,1);
	% X = [x1,y1,dotx1,doty1]
	x1    = X(1); % x position
	y1    = X(2); % y position
	dotx1 = X(3); % x velocity
	doty1 = X(4); % y velocity
	tmp1  = (dotx1^2+doty1^2);
	ddotx1 = (1/L^2)*(-g*x1*y1-x1*tmp1);
	ddoty1 = (1/L^2)*(g*x1^2-y1*tmp1);	
	dX(1) = dotx1;
	dX(2) =  doty1;
	dX(3) = ddotx1;
	dX(4) = ddoty1;
end