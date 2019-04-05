% Event for Halo orbit determinition

function [position,isterminal,direction] = EventHalo(t,q)
position = q(2); % The value to be zero  --> q(2) is y: x-z plane cross
isterminal = 1;  % Halt integration 
direction = 0;   % The zero can be approached from either direction