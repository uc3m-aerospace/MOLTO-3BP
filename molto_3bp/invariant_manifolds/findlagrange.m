function eqn = findlagrange(x,mustar)
%This function contains the equation used to find the Lagrangian points.
 eqn = -((1-mustar)/abs(mustar+x)^3)*(mustar+x)+(mustar/abs(1-mustar-x)^3)*(1-mustar-x) +x;
end
