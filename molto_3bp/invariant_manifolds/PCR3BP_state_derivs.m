function [derivs] = PCR3BP_state_derivs (et,state, params)
% This function computes the acceleration acting at time "et" on
% the fourth body according to the PBRFBP model.
% INPUTS:
%  - et     : Current non-dimensional time [sec]
%  - state  : Current non-dimensional S/C state in the synodic reference frame
%             [-,-]
% OUTPUTS:
%  - derivs: Total derivatives of the S/C state vector
%     - derivs (1:2) : Position derivatives
%     - derivs (3:4) : Velocity derivatives
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize derivative vector (4 components)
derivs = zeros(4,1);

% Get state
x = state(1);
y = state(2);
xdot = state(3);
ydot = state(4);

% Get parameters
mu1  = params.mu1;
mu2   = params.mu2;

% Intermediate variables
% r1 = sqrt((x+mu2)^2 + y^2);
% r2 = r2 = sqrt((x-mu1)^2 + y^2);
% Ubar = -1/2 * (mu1*r1^2 + mu2*r2^2) - mu1/r1 - mu2/r2
% dUbardx = simplify(diff(Ubar,x))
% dUbardy = simplify(diff(Ubar,y))

dUbardx = mu2*(mu1 - x) - mu1*(mu2 + x) - (mu2*(mu1 - x))/((mu1 - x)^2 + y^2)^(3/2) + (mu1*(mu2 + x))/((mu2 + x)^2 + y^2)^(3/2);
dUbardy = (mu1*y)/((mu2 + x)^2 + y^2)^(3/2) - mu2*y - mu1*y + (mu2*y)/((mu1 - x)^2 + y^2)^(3/2);

% Position derivatives
derivs(1:2) = state(3:4);

% Velocity derivatives
derivs(3) = 2*state(4)-dUbardx;
derivs(4) = -2*state(3)-dUbardy;

end
