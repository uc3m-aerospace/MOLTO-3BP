function [lookfor stop direction] = earth_EM_crossing_detection(et,states,params)

% The next lines of code define what is the variable that reaches 0 at the stop condition
x_EM = states(1);
y_EM = states(2);

% Moon position in EM reference system
xM_EM = 1-params.mu_EM;

if abs(x_EM-xM_EM) < 0.1 && abs(y_EM)<1
        variable_to_be_checked = x_EM-xM_EM;
else
    
    % State in SE reference system
    states_SE = EM2SE(states, et, params.mu_EM, params.mu_SE, params.phi_EM_0, params.L_EM, params.L_SE, params.T_EM, params.T_SE);
    x_SE = states_SE(1);
    y_SE = states_SE(2);
    
    % Earth position in SE reference system
    xE_SE = 1-params.mu_SE;
    
    % Earth position in EM reference system
    xE_EM = -params.mu_EM;
    
    states_Moon = EM2SE([xM_EM;0;0;0],et, params.mu_EM, params.mu_SE, params.phi_EM_0, params.L_EM, params.L_SE, params.T_EM, params.T_SE);
    xM_SE = states_Moon(1);
    yM_SE = states_Moon(2);
    
    tol = 10^-6;
    
    variable_to_be_checked = double(abs(x_SE-xE_SE) > tol) + double(y_SE < 3.5e-3);
    
end

% % Standard Matlab variable name that defines the stop condition (lookfor = 0)
lookfor = variable_to_be_checked;
% Next variable (stop) can be always set to 1
stop    = 1;                          % It means "stop" when event is located
% Next variable defines the direction of 0 crossing that determines a propagation stop
direction = 0;                        % The "lookfor" variable derivative at stop
%  1 : positive derivative
% -1 : negative derivative
%  0 : whatever derivative
end