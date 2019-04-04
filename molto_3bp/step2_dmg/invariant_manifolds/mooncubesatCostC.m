function [Mayer,Lagrange]=mooncubesatCostC(solcost,setup)

t0 = solcost.initial.time;
x0 = solcost.initial.state;
tf = solcost.terminal.time;
xf = solcost.terminal.state;
t  = solcost.time;
x  = solcost.state;
u  = solcost.control;
p  = solcost.parameter;
m = solcost.state(:,3);
TT = u(:,1)*setup.CONSTANTS.Tmax;

Mayer = 0;%-xf(3);
Lagrange = 0.5*(TT./m).^2;%zeros(size(t));