function events = mooncubesateventfun(solevents,setup)

xf = solevents.terminal.state;
x0 = solevents.initial.state;

mu_EM = setup.CONSTANTS.mu;

rf = sqrt((xf(1)-1+mu_EM)^2 + xf(2)^2);
Vm = sqrt(mu_EM/rf);
ang = atan2(xf(2),xf(1)-1+mu_EM);

ef = sqrt(((xf(4) - xf(2))+Vm*sin(ang))^2 + ((xf(5) + xf(1) + mu_EM -1)-Vm*cos(ang))^2);

% vf = sqrt((xf(4) - xf(2))^2 + (xf(5) + xf(1) + mu_EM -1)^2);
% ve = sqrt(2*mu_EM/rf);
% 
% ef = vf-ve;

e0_dir = (x0(1)+mu_EM)*(x0(4)-x0(2)) + x0(2)*(x0(5)+x0(1)+mu_EM);

rf = sqrt((xf(1)-1+mu_EM)^2 + xf(2)^2);

events=[ ef; e0_dir; rf];

