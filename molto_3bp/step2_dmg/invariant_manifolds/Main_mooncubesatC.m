%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%           OPTIMIZATION PROBLEM       %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Mooncubesatdae contains de differential equations

%%Mooncubesatfun contains the events

    %%% 0. Initialize General variables
           %a)Fixed parameters 
           %b)Smallsat parameters
           %c)Departure orbit parameters
           %d)Load initial guess
    %%% 1. GPOPS constraints
           %1.1 Boundary conditions 
           %1.2 Path constraints   
    %%% 2. GPOPS solver       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%% 0. Initialize General variables
clear all;
close all;
clc;
% a) Fixed Parameters
global CONSTANTS
G = 6.67408*10^-20;                             % Univ. grav cte[km3/(kgs2)]
mu_earth = 398600.435436096;
mu_sun = 1.32712440041939E+11;
mu_moon = 4902.80006616380;
m_earth = mu_earth/G;
m_sun = mu_sun/G;
m_moon = mu_moon/G;
l_em = 384405000;
mu_EM = mu_moon/(mu_moon+mu_earth);
mu1_EM = 1-mu_EM;
mu2_EM = mu_EM;
L_EM = 3.85*10^5;
T_EM = 2.361*10^6;
DU    = l_em;                                        % [m] Distance unit
TU    = 4.34811305*24*3600;                          % [days] Time unit
VU    = 1.02323281*10^3;                             % [m/s] Speed unit
MU    = mu_earth+mu_moon;                            % [kg]  Mass unit
CONSTANTS.rho = 3.88811143*10^2;               % Scaled Sun-(Earth+Moon) distance
CONSTANTS.mu = mu_moon/(mu_earth+mu_moon);     % Earth-Moon mass parameter [-]
CONSTANTS.m_s = mu_sun/(mu_earth+mu_moon);     % Scaled mass of the Sun
CONSTANTS.om_s = -9.25195985*10^(-1);          % Scaled angular velocity of the Sun
CONSTANTS.g0 = 9.81/(DU/TU^2);                 % Scale gravity acceleration on Earth surface
R_m   = 1738;                                  % [km] Mean Moon's radius

% c)Smallsat parameters
CONSTANTS.Tmax = 0.0003/(MU*DU/TU^2);           % Scaled maximum thrust
CONSTANTS.Isp = 2000/TU;                        % Scaled specific impulse
m00 = 2;                                        %Initial mass [kg]
% h_m   = 1000;        %[km] Altitude of arrival orbi

% c)Departure orbit parameters
CONSTANTS.phi_EM_0 =  225/180*pi;
R_GEO = 42164000/DU;
V_GEO = sqrt((1-CONSTANTS.mu)/R_GEO);

params.rho = CONSTANTS.rho;
params.mu = CONSTANTS.mu;
params.m_s = CONSTANTS.m_s;
params.om_s = CONSTANTS.om_s;
params.phi_EM_0 = CONSTANTS.phi_EM_0;

% d) Load initial guess:
load('.\Initial_guess\states_EM_complete', 'states_EM_complete');
load('.\Initial_guess\times_EM_complete', 'times_EM_complete');

% Change scaling of trajectory
IG_states(1,:) = states_EM_complete(1,1:round(end-1));
IG_states(2,:) = states_EM_complete(2,1:round(end-1));
IG_states(3,:) = states_EM_complete(3,1:round(end-1))*TU/(T_EM/(2*pi));
IG_states(4,:) = states_EM_complete(4,1:round(end-1))*TU/(T_EM/(2*pi));

IG_times = times_EM_complete(1:round(end-1))*(T_EM/(2*pi))/TU; 

npoints = 900;
nodes = 900;

% Resample IG times and states
idx = round(linspace(1,size(IG_times,2),npoints));
IG_times_n = IG_times(1,idx);
IG_states_n = IG_states(:,idx);

S0 = IG_states_n(:,1);
SF = IG_states_n(:,end); 
%% 1. GPOPS constraints 

%1.1 Boundary conditions
% a) Initial and final states obtained from initial guess (scaled)
x0 = S0(1);
xf = SF(1);
y0 = S0(2);
yf = SF(2);
xdot0 = S0(3);
ydot0 = S0(4);
% xdotf = SF(3);
% ydotf = SF(4);
m0    = m00/MU;
mf    = m00/MU;

% b)Change scaling initial and final states
r0 = sqrt((x0+CONSTANTS.mu)^2 + y0^2);
y0 = y0/r0*(R_GEO);
x0 = (x0+CONSTANTS.mu)/r0*R_GEO - CONSTANTS.mu;
r0 = sqrt((x0+CONSTANTS.mu)^2 + y0^2);
v0mag = sqrt((S0(3)-S0(2))^2+ (S0(4)+S0(1)+CONSTANTS.mu)^2);
ang0 = atan2(S0(2),S0(1)+CONSTANTS.mu);
xdot0 = -v0mag*sin(ang0)+y0;
ydot0 = v0mag*cos(ang0)-x0-CONSTANTS.mu;
rf = sqrt((xf-1+CONSTANTS.mu)^2 + yf^2);
yf = yf/rf*(3000+R_m)*1000/DU;
xf = (xf-1+CONSTANTS.mu)/rf*(3000+R_m)*1000/DU + 1 -CONSTANTS.mu;
rf = sqrt((xf-1+CONSTANTS.mu)^2 + yf^2);
ap = (rf + (500+R_m)*1000/DU)/2;
en = -CONSTANTS.mu/(2*ap);
vf = sqrt((en+CONSTANTS.mu/rf)*2);
ang = atan2(yf,(xf-1+CONSTANTS.mu));
xdotf = -vf*sin(ang)+yf;
ydotf = vf*cos(ang)-xf+1-CONSTANTS.mu;

IG_states_n(1,end) = xf;
IG_states_n(2,end) = yf;
IG_states_n(3,end) = xdotf;
IG_states_n(4,end) = ydotf;
IG_states_n(1,1) = x0;
IG_states_n(2,1) = y0;
IG_states_n(3,1) = xdot0;
IG_states_n(4,1) = ydot0;

% IG_states_n = IG_states_n(:,1:end-4);
% IG_times_n = IG_times_n(1:end-4);
% IG_states_n(:,end+1) = [xf;yf;xdotf;ydotf];
% IG_times_n(end+1) = IG_times(end);%+(IG_times_n(end)-IG_times_n(end-1));

% c) Initial boundary conditions
x0min = -R_GEO;
x0max = R_GEO;
y0min = -R_GEO;
y0max = R_GEO;
xdot0min = -4;
xdot0max = 4;
ydot0min = -4;
ydot0max = 4;
m0min = m00/MU;
m0max = m00/MU;

% d) Final boundary conditions
mfmin = 0;
mfmax = m0;
xfmin = (1-mu_EM)-10e6/DU; 
xfmax = (1-mu_EM)+10e6/DU;
yfmin = -10e6/DU;
yfmax = 10e6/DU;
xdotfmin = -4;
xdotfmax = 4;
ydotfmin = -4;
ydotfmax = 4;

% e) Time boundary conditions
t0 = IG_times_n(1);
t0min = t0;
t0max = t0;
tf = IG_times_n(end);
if tf<0
    tfmin = 2*IG_times_n(end);
    tfmax = 0.3*IG_times_n(end);
else
    tfmin = 0.3*IG_times_n(end);
    tfmax = 2*IG_times_n(end);
end

%% 1.2 Path constraints 

% a) States bounding
xmin = -5;
xmax = 5;
ymin = -5;
ymax = 5;
xdotmin = -4;
xdotmax = 4;
ydotmin = -4;
ydotmax = 4;
mmin = 0;
mmax = m0;

% b) Control variables bounding
throttlemin =   0;
throttlemax =   1;
alphamin = 0;
alphamax = 2*pi;

%% GPOPS solver
% 2.1 Phase 1 Information
iphase = 1;
%Limits
limits(iphase).nodes          = nodes;
limits(iphase).time.min       = [t0min tfmin];
limits(iphase).time.max       = [t0max tfmax];
limits(iphase).state.min(1,:) = [x0min xmin xfmin];
limits(iphase).state.max(1,:) = [x0max xmax xfmax];
limits(iphase).state.min(2,:) = [y0min ymin yfmin];
limits(iphase).state.max(2,:) = [y0max ymax yfmax];
limits(iphase).state.min(3,:) = [m0 mmin mfmin];
limits(iphase).state.max(3,:) = [m0 mmax mfmax];
limits(iphase).state.min(4,:) = [xdot0min xdotmin xdotfmin];
limits(iphase).state.max(4,:) = [xdot0max xdotmax xdotfmax];
limits(iphase).state.min(5,:) = [ydot0min ydotmin ydotfmin];
limits(iphase).state.max(5,:) = [ydot0max ydotmax ydotfmax];

limits(iphase).control.min(1,:)    = throttlemin;
limits(iphase).control.max(1,:)    = throttlemax;
limits(iphase).control.min(2,:)    = alphamin;
limits(iphase).control.max(2,:)    = alphamax;

limits(iphase).parameter.min  = [];
limits(iphase).parameter.max  = [];

limits(iphase).path.min       = [];
limits(iphase).path.max       = [];

limits(iphase).duration.min   = [];
limits(iphase).duration.max   = [];

%Define upper and lower limits of the events defined in mooncubesateventfun
limits(iphase).event.min   = [-1e-20;(100 + R_m)*1000/DU;-1e-20;0*sqrt((1-CONSTANTS.mu)/r0);R_GEO+1e-20];
limits(iphase).event.max   = [1e-20; (8000 + R_m)*1000/DU;1e-20;0.5*sqrt((1-CONSTANTS.mu)/r0);R_GEO-1e-20];

%Initial Guess
guess(iphase).time            = IG_times_n';
guess(iphase).state(:,1)      = IG_states_n(1,:);
guess(iphase).state(:,2)      = IG_states_n(2,:);
guess(iphase).state(:,3)      = m0*ones(size(IG_states_n(1,:)));
guess(iphase).state(:,4)      = IG_states_n(3,:);
guess(iphase).state(:,5)      = IG_states_n(4,:);
guess(iphase).control(:,1)    = zeros(size(IG_states_n(1,:)));
guess(iphase).control(:,2)    = zeros(size(IG_states_n(1,:)));
guess(iphase).parameter       = [];

linkages = [];
setup.CONSTANTS = CONSTANTS;
setup.name  = 'mooncubesatProblem';
setup.funcs.cost = 'mooncubesatCostC';
setup.funcs.dae = 'mooncubesatDaeC';
setup.funcs.event = 'mooncubesateventfunC';
setup.limits = limits;
setup.guess = guess;
setup.linkages = linkages;
setup.derivatives = 'numerical'; %numerical/complex/analytical
setup.direction = 'increasing';
setup.parallel = 'no';
setup.autoscale = 'off'; %on/off
setup.solver ='ipopt';% snopt /ipopt
setup.method ='collocation'; %pseudospectral/collocation

% 2.2 Run & save the solution
output = DMG(setup);

solution = output.solution;
t = now;
date = datestr(t,'mmmm_dd_HH_MM_SS');
save(strcat([date '_solution']),'solution')

%% 3. Postprocessing

%Thrust
figure
plot(solution.time*TU/3600/24,solution.control(:,1)*CONSTANTS.Tmax*(MU*DU/TU^2),'b','Linewidth',2)
xlabel('Time [days]','interpreter','latex','fontsize',18);
ylabel('Thrust [N]','interpreter','latex','fontsize',18);
grid on
set(gca, 'ticklabelinterpreter', 'latex', 'fontsize', 15)
savefig(strcat([date 'thrust.fig']))

%Thrust direction
figure
plot(solution.time*TU/3600/24,solution.control(:,2)*180/pi,'b','Linewidth',2)
xlabel('Time [days]','interpreter','latex','fontsize',18);
ylabel('Angle [deg]','interpreter','latex','fontsize',18);
grid on
set(gca, 'ticklabelinterpreter', 'latex', 'fontsize', 15)
savefig(strcat([date 'angle.fig']))

%Trajectory
figure
plot(solution.state(:,1),solution.state(:,2),'b','Linewidth',2)
hold on
plot(IG_states_n(1,:),IG_states_n(2,:),'r','Linewidth',2)
xlabel('$X_{EM}$','interpreter','latex','fontsize',18);
ylabel('$Y_{EM}$','interpreter','latex','fontsize',18);
grid on
plot(-CONSTANTS.mu,0,'ko','markerfacecolor','k','markersize',4)%interp2(X,Y,U,-mu2,0),'bo')
 text(-CONSTANTS.mu-0.17,0,'$m_E$','interpreter','latex','fontsize',15)
plot(1-CONSTANTS.mu,0,'ko','markerfacecolor','k','markersize',4)%interp2(X,Y,U,-mu2,0),'ro')
 text(1-CONSTANTS.mu-0.17,0,'$m_M$','interpreter','latex','fontsize',15)
set(gca, 'ticklabelinterpreter', 'latex', 'fontsize', 15)
savefig(strcat([date 'trajectory.fig']))

%Mass
figure
plot(solution.time*TU/3600/24,solution.state(:,3)*MU,'b','Linewidth',2)
hold on
xlabel('Time [days]','interpreter','latex','fontsize',18);
ylabel('Mass [kg]','interpreter','latex','fontsize',18);
grid on
set(gca, 'ticklabelinterpreter', 'latex', 'fontsize', 15)
savefig(strcat([date 'mass.fig']))

params.rho = CONSTANTS.rho;
params.mu = CONSTANTS.mu;
params.m_s = CONSTANTS.m_s;
params.om_s = CONSTANTS.om_s;
params_EM.phi_EM_0 = 0;

S0 = solution.state(end,[1:2 4:5])';
deltat = 0.5;
et0 = solution.time(end);
prnt_out_dt = 0.00005;
stop_fun = 'none';
[SF, etf, states, times] = PCRFBP_propagator (S0, params,...
    et0, deltat, prnt_out_dt, stop_fun);


figure
plot(states(1,:),states(2,:),'b','linewidth',2)
hold on
plot(1-CONSTANTS.mu,0,'ko','markerfacecolor','k','markersize',4)%interp2(X,Y,U,-mu2,0),'ro')
 text(1-CONSTANTS.mu-0.002,0,'$m_M$','interpreter','latex','fontsize',15)
grid on
xlabel('$X_{EM}$','interpreter','latex','fontsize',18);
ylabel('$Y_{EM}$','interpreter','latex','fontsize',18);
plot(states(1,1),states(2,1),'bo','markerfacecolor','b','markersize',6)
%set(gca,'ticklabelinterpreter','latex','fontsize',15)
%savefig(strcat([date 'moon.fig']))

x0sol = solution.state(1,1);
y0sol = solution.state(1,2);
xdot0sol = solution.state(1,4);
ydot0sol = solution.state(1,5);

r0sol = sqrt((x0sol+CONSTANTS.mu)^2 + y0sol^2);
ang0 =  atan2(y0sol,x0sol+mu_EM);
V_c0 = sqrt((1-CONSTANTS.mu)/r0sol);
dv0 = sqrt((xdot0sol-y0sol+V_c0*sin(ang0))^2 + (ydot0sol+x0sol+CONSTANTS.mu-V_c0*cos(ang0))^2);

r_moon = sqrt((states(1,:)-1+CONSTANTS.mu).^2 + states(2,:).^2);

figure
plot(times*TU/3600/24,r_moon*DU/1000-R_m,'b','Linewidth',2)
grid on
xlabel('Time [s]','interpreter','latex','fontsize',18);
ylabel('$h_{Moon}$ [km]','interpreter','latex','fontsize',18);

figure
X_ECI = (solution.state(:,1) + CONSTANTS.mu).*cos(solution.time') - solution.state(:,2).*sin(solution.time');
Y_ECI = (solution.state(:,1) + CONSTANTS.mu).*sin(solution.time') + solution.state(:,2).*cos(solution.time');
VX_ECI = (solution.state(:,4) - solution.state(:,2)).*cos(solution.time') - (solution.state(:,5) + solution.state(:,1) + CONSTANTS.mu).*sin(solution.time');
VY_ECI = (solution.state(:,4) - solution.state(:,2)).*sin(solution.time') - (solution.state(:,5) + solution.state(:,1) + CONSTANTS.mu).*cos(solution.time');

plot(X_ECI,Y_ECI)
hold on
quiver(X_ECI,Y_ECI,VX_ECI,VY_ECI)