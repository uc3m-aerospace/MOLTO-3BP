function [f,cons] = fitness_nsga2(x,params)

%% Identify decision variables
alpha = x(1);
beta = x(2);
ti = x(3);
delta = x(4);
%
r0 = params.r0;
rf = params.rf;
mu = params.mu;
DU = params.DU;
VU = params.VU;
%
prnt_out_dt = 0.001;
stop_fun = 'none';

%% Initial state
x0    = r0*cos(alpha)-mu;
y0    = r0*sin(alpha);
v0    = beta*sqrt((1-mu)/r0);
x0dot = -(v0-r0)*sin(alpha);
y0dot = (v0-r0)*cos(alpha);

S0 = [x0;y0;x0dot;y0dot];

%% Propagate PCRFBP
[SF, etf, states, times] = PCRFBP_propagator (S0, params,...
    ti, delta, prnt_out_dt, stop_fun);

xf = SF(1);
yf = SF(2);
xfdot = SF(3);
yfdot = SF(4);
alpha = atan2(yf,xf+mu-1);
Vm = sqrt(mu/rf);

%% Evaluate objective functions

f(1) = abs(sqrt((xf+mu-1)^2+yf^2)-rf)*DU/1000; %delta r (km - dimensional units)
f(2) = abs(sqrt((-Vm*sin(alpha)-(xfdot-yf))^2 + (Vm*cos(alpha)-(yfdot+xf+mu-1))^2))*VU/1000; % velocity vector error (km/s - dimensional units=

cons = 0;

