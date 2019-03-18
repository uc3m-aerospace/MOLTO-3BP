function [SF, etf, states, times] = PCRFBP_propagator (S0, params,...
                                    et0, deltat, prnt_out_dt, stop_fun)
% ------------------- NUMERICAL INTEGRATION FUNCTION ----------------------
% This function propagates the state of a S/C according to the PBRFBP.
%
% -------------------------------------------------------------------------

% Initial S/C state vector in the input inertial reference frame
state = S0;
% Final ephemeris time
etf = et0 + deltat;

% Function that computes the state derivatives
derivs = @(et,state)PCRFBP_state_derivs (et, state, params);

% % Options of the Runge Kutta solver
% % Maximum integration time step is 1/50 of the orbital period
% sc_elems = cspice_oscelt( state(1:6), et0, othr_params.mu_cen );
% sma      = sc_elems(1) / (1-sc_elems(2));
% if (sma > 0)
%     orb_per  = 2*pi*sqrt(sma^3/othr_params.mu_cen);
%     max_step = 1/10 * orb_per;
% else
    max_step = 1e20; % 1/50 Moon orbit of 100 km altitude
% end

% Additional exit condition of the Cowell's propagator (if it is reached
% before the final propagation time)
options = odeset('RelTol',1e-9,'AbsTol',1e-9, 'MaxStep',max_step);
% Set events function to the input function handle
if ( strcmp(stop_fun, 'none') == 0 )
    options.Events = stop_fun;
end

% Avoid errors when the print out step has been set higher than the
% propagation duration
if (prnt_out_dt > etf-et0)
    prnt_out_dt = etf-et0;
end
% ---------- SOLVE FOR THE TRAJECTORY WITH AN ODE45 INTEGRATOR ------------
[times,states] = ode45(derivs,et0:prnt_out_dt:etf,state,options); 

% When propagating in time, assure that the final propagation time is 
% included in the simulation
if (strcmp(stop_fun, 'none') == 1 && times(end) < etf)
    time0 = times(end);
    timef = etf;
    dt = timef-time0;
    S0 = states(end,1:4)';
    [ttt,sss] = ode45(derivs,time0:dt:timef,S0,options);
    times = [times;ttt(end)];
    states = [states;sss(end,1:4)];
end

% Reformulate output vector dimensions
times  = times';
states = states';

% Update the final S/C state value and ephemeris time
SF = states(1:4,end);
etf = times(end);

return;
