function daeout = mooncubesatDae(soldae,setup)

t = soldae.time;
state = soldae.state;
u = soldae.control;
   
throttle = u(:,1);
alpha = u(:,2);

% Get parameters
rho  = setup.CONSTANTS.rho;
mu   = setup.CONSTANTS.mu;
m_s  = setup.CONSTANTS.m_s;
om_s = setup.CONSTANTS.om_s;
Isp  = setup.CONSTANTS.Isp;
g0   = setup.CONSTANTS.g0;
Tmax = setup.CONSTANTS.Tmax;

% Get state
x = state(:,1);
y = state(:,2);
m = state(:,3);

% Intermediate variables
% r1 = sqrt((x+mu)^2+y^2); % Distance SC from Earth
% r2 = sqrt((x+mu-1)^2+y^2); % Distance SC from Moon
% r3 = sqrt((x-rho*cos(om_s*t))^2 + (y-rho*sin(om_s*t))^2); % Distance SC from Sun
% Om3 = 0.5*(x^2+y^2) + (1-mu)/r1 + mu/r2 + 0.5*mu*(1-mu);
% Om4 = Om3 + m_s/r3 - m_s/rho^2 * (x*cos(om_s*t)+y*sin(om_s*t));

dOm4dx = x - (m_s.*cos(om_s.*t))./rho.^2 - (mu.*(mu + x - 1))./((mu + x - 1).^2 + y.^2).^(3/2) + ((2.*mu + 2.*x).*(mu - 1))./(2.*((mu + x).^2 + y.^2).^(3/2)) - (m_s.*(2.*x - 2.*rho.*cos(om_s.*t)))./(2.*((x - rho.*cos(om_s.*t)).^2 + (y - rho.*sin(om_s.*t)).^2).^(3/2));
dOm4dy = y - (m_s.*sin(om_s.*t))./rho.^2 - (mu.*y)./((mu + x - 1).^2 + y.^2).^(3/2) - (m_s.*(2.*y - 2.*rho.*sin(om_s.*t)))./(2.*((x - rho.*cos(om_s.*t)).^2 + (y - rho.*sin(om_s.*t)).^2).^(3/2)) + (y.*(mu - 1))./((mu + x).^2 + y.^2).^(3/2);

% Position derivatives
xdot = state(:,4);
ydot = state(:,5);

% Mass derivative
T = throttle * Tmax;
mdot = -T/(Isp*g0);

% Velocity derivatives
Tx = T.*cos(alpha);
Ty = T.*sin(alpha);
xdotdot = dOm4dx + 2*ydot + Tx./m;
ydotdot = dOm4dy - 2*xdot + Ty./m;

daeout = [xdot ydot mdot xdotdot ydotdot];

