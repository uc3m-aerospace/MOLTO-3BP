% Compute derivatives of Om4

syms x y mu m_s rho t om_s

r1 = sqrt((x+mu)^2+y^2); % Distance SC from Earth
r2 = sqrt((x+mu-1)^2+y^2); % Distance SC from Moon
r3 = sqrt((x-rho*cos(om_s*t))^2 + (y-rho*sin(om_s*t))^2); % Distance SC from Sun

Om3 = 0.5*(x^2+y^2) + (1-mu)/r1 + mu/r2 + 0.5*mu*(1-mu);

Om4 = Om3 + m_s/r3 - m_s/rho^2 * (x*cos(om_s*t)+y*sin(om_s*t));

dOm4dx = diff(Om4,x);
dOm4dy = diff(Om4,y);