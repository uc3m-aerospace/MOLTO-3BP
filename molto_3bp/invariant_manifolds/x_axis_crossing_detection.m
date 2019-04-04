function [lookfor stop direction] = x_axis_crossing_detection(et,states_aux)

% The next lines of code define what is the variable that reaches 0 at the stop condition
x = states_aux(1);
y = states_aux(2);
variable_to_be_checked = y;
% % Standard Matlab variable name that defines the stop condition (lookfor = 0)
lookfor = variable_to_be_checked;
% Next variable (stop) can be always set to 1
stop    = 1;                          % It means "stop" when event is located
% Next variable defines the direction of 0 crossing that determines a propagation stop
direction = 0;                        % The "lookfor" variable derivative at stop
                                      %  1 : positive derivative
									  % -1 : negative derivative
									  %  0 : whatever derivative