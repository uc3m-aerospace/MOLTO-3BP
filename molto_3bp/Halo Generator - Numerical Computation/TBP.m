function dqdt = TBP(t,q)
% Control -> Choice case 

opt=2; % 1 -> Sun - Earth+Moon System
       % 2 -> Earth - Moon Syst

%% Init

x=q(1);
y=q(2);
z=q(3);
vx=q(4);
vy=q(5);
vz=q(6);

%% Data Required
%INITIAL DATA

mE=5.9722e24;
mM=7.342e22;
mS=1.98845114e30;

if opt==1    % Sun-Earth
    mu=(mE+mM)/(mE+mM+mS);
else 
    mu=0.012150585609624;%1.215e-2;
end

mu1 = 1-mu;
mu2 = mu;
r1 = ((x+mu2)^2+y^2+z^2)^(1/2);
r2 = ((x-mu1)^2+y^2+z^2)^(1/2);

% U

Ux = (mu2*(2*mu1 - 2*x))/2 - (mu1*(2*mu2 + 2*x))/2 ...
    - (mu2*(2*mu1 - 2*x))/(2*((mu1 - x)^2 + y^2 + z^2)^(3/2)) ...
    + (mu1*(2*mu2 + 2*x))/(2*((mu2 + x)^2 + y^2 + z^2)^(3/2));
Uy = (mu1*y)/((mu2 + x)^2 + y^2 + z^2)^(3/2) - mu2*y - mu1*y ...
    + (mu2*y)/((mu1 - x)^2 + y^2 + z^2)^(3/2);
Uz = (mu1*z)/((mu2 + x)^2 + y^2 + z^2)^(3/2) + ...
    (mu2*z)/((mu1 - x)^2 + y^2 + z^2)^(3/2);

%% Ode Work
xdot = vx;
ydot = vy;
zdot = vz;
xdotdot =  2*vy-Ux;
ydotdot = -2*vx-Uy;
zdotdot = -Uz;

dqdt=zeros(size(q));
dqdt(1)=xdot;
dqdt(2)=ydot;
dqdt(3)=zdot;
dqdt(4)=xdotdot;
dqdt(5)=ydotdot;
dqdt(6)=zdotdot;
