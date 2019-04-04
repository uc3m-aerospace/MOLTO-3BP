clear all
close all
clc

% Load initial guess
path = './';
addpath(path);
name = 'solution_April_13_22_59_57_states.txt';
filename = strcat([path name]);
delimiterIn = ' ';
IG_states = importdata(filename, delimiterIn);
S0 = IG_states(1:4,1);
SF = IG_states(1:4,end);
IG_times = IG_states(5,:);
IG_states = IG_states(1:4,:);

npoints=100;

%% Resample IG times and states
idx = round(linspace(1,size(IG_times,2),npoints));
IG_times_n = IG_times(1,idx);
IG_states_n = IG_states(:,idx);

global CONSTANTS

mu_earth = 398600.435436096;
mu_sun = 1.32712440041939E+11;
mu_moon = 4902.80006616380;
l_em  = 384405000;
DU    = l_em;                                  % [m] Distance unit
TU    = 4.34811305*24*3600;                    % [days] Time unit
VU    = 1.02323281*10^3;                       % [m/s] Speed unit
MU    = mu_earth+mu_moon;                      % [kg]  Mass unit
CONSTANTS.rho = 3.88811143*10^2;               % Scaled Sun-(Earth+Moon) distance
CONSTANTS.mu = mu_moon/(mu_earth+mu_moon);     % Earth-Moon mass parameter [-]
CONSTANTS.m_s = mu_sun/(mu_earth+mu_moon);     % Scaled mass of the Sun
CONSTANTS.om_s = -9.25195985*10^(-1);          % Scaled angular velocity of the Sun
CONSTANTS.Tmax = 0.0003/(MU*DU/TU^2);          % Scaled maximum thrust
CONSTANTS.Isp = 2000/TU;                       % Scaled specific impulse
CONSTANTS.g0 = 9.81/(DU/TU^2);                 % Scale gravity acceleration on Earth surface
R_m   = 1738;                                  % [km] Mean Moon's radius

x0 = S0(1);
xf = SF(1);
y0 = S0(2);
yf = SF(2);
xdot0 = S0(3);
xdotf = SF(3);
ydot0 = S0(4);
ydotf = SF(4);
m0    = 2/MU;
mf    = 2/MU;
m0min = 2/MU;
m0max = 2/MU;
mfmin = 0;
mfmax = m0;

rf = sqrt((xf-1+CONSTANTS.mu)^2 + yf^2);
vf = sqrt(CONSTANTS.mu/rf);
ang = atan2(yf,(xf-1+CONSTANTS.mu));
xdotf = -vf*sin(ang)+yf;
ydotf = vf*cos(ang)-xf+1-CONSTANTS.mu;

xdot0min = min([xdot0,0.99*xdot0,1.01*xdot0]);
xdot0max = max([xdot0,0.99*xdot0,1.01*xdot0]);
ydot0min = min([ydot0,0.99*ydot0,1.01*ydot0]);
ydot0max = max([ydot0,0.99*ydot0,1.01*ydot0]);
xfmin = (1-CONSTANTS.mu)-3e6/DU;
xfmax = (1-CONSTANTS.mu)+3e6/DU;
yfmin = -3e6/DU;
yfmax = 3e6/DU;
xdotfmin = -2;
xdotfmax = 2;
ydotfmin = -2;
ydotfmax = 2;

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

throttlemin =   0;
throttlemax =   1;
alphamin = 0;
alphamax = 2*pi;

t0 = IG_times_n(1);
t0min = t0;
t0max = t0;
tf = IG_times_n(end);
tfmin = 0.8*IG_times(end);
tfmax = 1.2*IG_times(end);

% Phase 1 Information
iphase = 1;
limits(iphase).nodes           = npoints;
limits(iphase).time.min        = [t0min tfmin];
limits(iphase).time.max        = [t0max tfmax];
limits(iphase).state.min(1,:) = [x0 xmin xfmin];
limits(iphase).state.max(1,:) = [x0 xmax xfmax];
limits(iphase).state.min(2,:) = [y0 ymin yfmin];
limits(iphase).state.max(2,:) = [y0 ymax yfmax];
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
limits(iphase).event.min   = [-1e-20;-1e-20;(4000+R_m)/DU];
limits(iphase).event.max   = [1e-20;1e-20;(15000+R_m)/DU];

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
setup.funcs.cost = 'mooncubesatCost';
setup.funcs.dae = 'mooncubesatDae';
setup.funcs.event = 'mooncubesateventfun';
setup.limits = limits;
setup.guess = guess;
setup.linkages = linkages;
setup.derivatives = 'numerical'; %numerical/complex/analytical
setup.direction = 'increasing';
setup.parallel = 'no';
setup.autoscale = 'off'; %on/off
setup.solver ='ipopt';% snopt /ipopt
setup.method ='collocation'; %pseudospectral/collocation

output = DMG(setup);

solution = output.solution;
t = now;
date = datestr(t,'mmmm_dd_HH_MM_SS');
save(strcat([date '_solution']),'solution')

%% 
figure(1)
plot(solution.time*TU/3600/24,solution.control(:,1)*CONSTANTS.Tmax*(MU*DU/TU^2),'b','Linewidth',2)
xlabel('Time [days]','interpreter','latex','fontsize',18);
ylabel('Thrust [N]','interpreter','latex','fontsize',18);
grid on
%set(gca, 'ticklabelinterpreter', 'latex', 'fontsize', 15)
savefig(strcat([date 'thrust.fig']))

figure(2)
plot(solution.time*TU/3600/24,solution.control(:,2)*180/pi,'b','Linewidth',2)
xlabel('Time [days]','interpreter','latex','fontsize',18);
ylabel('$\ [deg]','interpreter','latex','fontsize',18);
grid on
%set(gca, 'ticklabelinterpreter', 'latex', 'fontsize', 15)
savefig(strcat([date 'angle.fig']))

figure(3)
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
%set(gca, 'ticklabelinterpreter', 'latex', 'fontsize', 15)
savefig(strcat([date 'trajectory.fig']))

figure(4)
plot(solution.time*TU/3600/24,solution.state(:,3)*MU,'b','Linewidth',2)
hold on
xlabel('Time [days]','interpreter','latex','fontsize',18);
ylabel('Mass [kg]','interpreter','latex','fontsize',18);
grid on
%set(gca, 'ticklabelinterpreter', 'latex', 'fontsize', 15)
savefig(strcat([date 'mass.fig']))

% figure(6)
% for ii = 1:size(solution.state,1)
%     states_EM_solution(1,:) = solution.state(:,1);
%     states_EM_solution(2,:) = solution.state(:,2);
%     states_EM_solution(3,:) = 
% 

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


figure(5)
plot(states(1,:),states(2,:),'b','linewidth',2)
hold on
plot(1-CONSTANTS.mu,0,'ko','markerfacecolor','k','markersize',4)%interp2(X,Y,U,-mu2,0),'ro')
text(1-CONSTANTS.mu-0.002,0,'$m_M$','interpreter','latex','fontsize',15)
grid on
xlabel('$X_{EM}$','interpreter','latex','fontsize',18);
ylabel('$Y_{EM}$','interpreter','latex','fontsize',18);
plot(states(1,1),states(2,1),'bo','markerfacecolor','b','markersize',6)
%set(gca,'ticklabelinterpreter','latex','fontsize',15)
savefig(strcat([date 'moon.fig']))

x0sol = solution.state(1,1);
y0sol = solution.state(1,2);
xdot0sol = solution.state(1,4);
ydot0sol = solution.state(1,5);


r0sol = sqrt((x0sol+CONSTANTS.mu)^2 + y0sol^2);
dv0 = sqrt((xdot0sol-y0sol)^2 + (ydot0sol+x0sol+CONSTANTS.mu)^2) - sqrt((1-CONSTANTS.mu)/r0sol);



