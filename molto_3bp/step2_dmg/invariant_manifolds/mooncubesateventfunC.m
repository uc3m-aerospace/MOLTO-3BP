function events = mooncubesateventfunC(solevents,setup)

xf = solevents.terminal.state;
x0 = solevents.initial.state;

mu_EM = setup.CONSTANTS.mu;

r0 = sqrt((x0(1)+mu_EM)^2 + x0(2)^2);
rf = sqrt((xf(1)-1+mu_EM)^2 + xf(2)^2);
Vm = sqrt(mu_EM/rf);
ang = atan2(xf(2),xf(1)-1+mu_EM);

ef1 = sqrt(((xf(4) - xf(2))+Vm*sin(ang))^2 + ((xf(5) + xf(1) + mu_EM -1)-Vm*cos(ang))^2);
ef2 = rf;

%ef2 = (x0(1) + mu_EM)*(x0(4) - x0(2)) + x0(2)*(x0(5)+x0(1)+mu);

% vf = sqrt((xf(4) - xf(2))^2 + (xf(5) + xf(1) + mu_EM -1)^2);
% ve = sqrt(2*mu_EM/rf);
% 
% ef = vf-ve;

% ef_dir = (xf(1)+mu_EM-1)*(xf(4)-xf(2)) + xf(2)*(xf(5)+xf(1)+mu_EM-1);

e0_dir = (x0(1)+mu_EM)*(x0(4)-x0(2)) + x0(2)*(x0(5)+x0(1)+mu_EM);

ang0 =  atan2(x0(2),x0(1)+mu_EM);
Vm0 = sqrt((1-mu_EM)/r0);

dv0_mag = sqrt((x0(4)-x0(2)+Vm0*sin(ang0))^2+(x0(5)+x0(1)+mu_EM-Vm0*cos(ang0))^2);

events=[ ef1; ef2; e0_dir; dv0_mag; r0];

