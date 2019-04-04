function [lookfor stop direction] = earth_SE_crossing_detection(et,states,params)

% The next lines of code define what is the variable that reaches 0 at the stop condition
x = states(1);
y = states(2);

tol = 10^-5;

if abs(x-params.mu1) < 10^-2 && (x-params.mu1) <0
    variable_to_be_checked = double(abs(y) > tol);
elseif abs(y) > 0.02
    variable_to_be_checked = 0;
else
    variable_to_be_checked = 1;
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