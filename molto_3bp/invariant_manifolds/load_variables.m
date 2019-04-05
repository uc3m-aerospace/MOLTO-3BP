%%Load parameters 
EM = load_spice_kernels(path2mice);

%%INITIAL DATA

% Gravitational constants:
mu_earth = cspice_bodvrd( 'EARTH', 'GM', 1 );   % Earth grav. pot [km3/s2]
mu_sun   = cspice_bodvrd( 'SUN', 'GM', 1 );     % Sun grav. pot [km3/s2]
mu_moon  = cspice_bodvrd( 'MOON', 'GM', 1);     % Moon grav. pot [km3/s2]
G = 6.67408*10^-20;       

%Earth-Moon system
mu_EM = mu_moon/(mu_moon+mu_earth);
mu1_EM = 1-mu_EM;
mu2_EM = mu_EM;

L_EM = 3.85*10^5;%distance between primaries (p 33)
T_EM = 2.361*10^6;%orbital period of primaries 27.32 days

%Sun-Earth system
mu_SE = mu_earth/(mu_earth+mu_sun);
mu1_SE = 1-mu_SE;
mu2_SE = mu_SE;

L_SE = 1.496*10^8; %Distance Sun-Earth
T_SE = 3.156*10^7; %Period of Earth around Sun: 365 days 
R_SE = 6378/L_SE;
r_Moon = 1738;
r_Moon_SE = r_Moon/L_SE;