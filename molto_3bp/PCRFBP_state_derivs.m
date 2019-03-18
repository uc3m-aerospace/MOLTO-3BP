function [derivs] = PCRFBP_state_derivs (t, state, params)
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
rho  = params.rho;
mu   = params.mu;
m_s  = params.m_s;
om_s = params.om_s;

% Intermediate variables
% r1 = sqrt((x+mu)^2+y^2); % Distance SC from Earth
% r2 = sqrt((x+mu-1)^2+y^2); % Distance SC from Moon
% r3 = sqrt((x-rho*cos(om_s*t))^2 + (y-rho*sin(om_s*t))^2); % Distance SC from Sun
% Om3 = 0.5*(x^2+y^2) + (1-mu)/r1 + mu/r2 + 0.5*mu*(1-mu);
% Om4 = Om3 + m_s/r3 - m_s/rho^2 * (x*cos(om_s*t)+y*sin(om_s*t));

dOm4dx = x - (m_s*cos(om_s*t))/rho^2 - (mu*(mu + x - 1))/((mu + x - 1)^2 + y^2)^(3/2) + ((2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(3/2)) - (m_s*(2*x - 2*rho*cos(om_s*t)))/(2*((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(3/2));
dOm4dy = y - (m_s*sin(om_s*t))/rho^2 - (mu*y)/((mu + x - 1)^2 + y^2)^(3/2) - (m_s*(2*y - 2*rho*sin(om_s*t)))/(2*((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(3/2)) + (y*(mu - 1))/((mu + x)^2 + y^2)^(3/2);

% Position derivatives
derivs(1:2) = state(3:4);

% Velocity derivatives
derivs(3) = dOm4dx + 2*state(4);
derivs(4) = dOm4dy - 2*state(3);

end
