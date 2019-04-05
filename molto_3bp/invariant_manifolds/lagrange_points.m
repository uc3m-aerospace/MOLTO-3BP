function pos = lagrange_points(mu,guess)
%
% Solve equation to obtain colinear Lagrange points
%
options = optimset('FunValCheck','on','TolX',10^-10,'TolFun',10^-10);
poslagrange3 = fzero(@(x) findlagrange(x,mu),guess(3),options);
poslagrange1 = fzero(@(x) findlagrange(x,mu),guess(1),options);
poslagrange2 = fzero(@(x) findlagrange(x,mu),guess(2),options);

r1 = -mu;
pos(1,:) = [poslagrange1,0];
pos(2,:) = [poslagrange2,0];
pos(3,:) = [poslagrange3,0];
pos(4,:) = [cos(pi/3)+r1,sin(pi/3)];
pos(5,:) = [cos(pi/3)+r1,-sin(pi/3)];

end