%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%      INITIAL GUESS (MANIFOLDS)   %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This code is divided in the following sections:

    %%% 0. Initialize General variables
    %%% 1. DEFINE STATES TO MATCH MANIFOLDS AT SECTION B
    %%% 2. EARTH-SCAPE LEG: Insert in SE model
    %%% 3. MOON-CAPTURE-LEG: Insert in EM model
    %%% 4. Complete trajectory     

%% 0. Initialize General variables
close all
clear all
clc 
%Mice package, Gravitational constants, Earth-Moon, and Sun-Earth constants
load_variables;

%Initial angle formed by x-axis of SE and EM system
phi_EM_0 = 225*pi/180;

phi_B = 90*pi/180;
%Define parameters (input variable to compute manifolds)
params_EM.mu1 = mu1_EM;
params_EM.mu2 = mu2_EM;
params_EM.mu_EM = mu2_EM;
params_EM.mu_SE = mu2_SE;
params_EM.L_EM = L_EM;
params_EM.L_SE = L_SE;
params_EM.T_EM = T_EM;
params_EM.T_SE = T_SE;
params_EM.phi_EM_0 = phi_EM_0;

params_SE.mu1 = mu1_SE;
params_SE.mu2 = mu2_SE;
%% 1. DEFINE STATES TO MATCH MANIFOLDS AT SECTION B

%1.1 Load states
load(strcat(['SF_s_EM_SE' num2str(phi_EM_0*180/pi)]),'SF_s_EM_SE')
load(strcat(['SF_u_SE_phiB_' num2str(phi_B*180/pi)]),'SF_u_SE_phiB')

% Final states of unstable manifolds at SB
x_u = SF_u_SE_phiB(1,:);
y_u = SF_u_SE_phiB(2,:);
xdot_u = SF_u_SE_phiB(3,:);
ydot_u = SF_u_SE_phiB(4,:);

%Stable manifolds EM converted into SE ref frame
x_s = SF_s_EM_SE(1,:);
y_s = SF_s_EM_SE(2,:);
xdot_s = SF_s_EM_SE(3,:);
ydot_s = SF_s_EM_SE(4,:);

% Match states at Section B: SE unstable with EM stable

idx_SE = find(y_u > 4.2e-3 & y_u < 4.24e-3)
x_u_idx = x_u(idx_SE);
y_u_idx = y_u(idx_SE);
xdot_u_idx = xdot_u(idx_SE);
ydot_u_idx = ydot_u(idx_SE);

%Match states at Section A: Search for the stable manifold states that start at 36000km from Earth (GEO)
idx_EM = find(y_s > 4.2e-3 & y_s < 4.24e-3);
x_s_idx = x_s(idx_EM);
y_s_idx = y_s(idx_EM);
xdot_s_idx = xdot_s(idx_EM);
ydot_s_idx = ydot_s(idx_EM);

%% 2. Insert in SE model

% Initial Conditions for integration => SO_SE
y0SE = [y_u_idx ;4.14e-3 ;4.145e-3; 4.15e-3;4.14e-3;4.16e-3;4.17e-3;4.18e-3;4.195e-3;4.197e-3;4.193e-3;4.194e-3;4.1937e-3];
ydot0SE = [(ydot_u_idx*9/10 + ydot_s_idx*1/10);-0.01435; -0.0143;-0.01425;-0.01435;-0.0142;-0.0141;-0.014;-0.0139;-0.0139;-0.0139;-0.0139;-0.0139];

x0_SE = mu1_SE;
y0_SE = y0SE(1)   %y_u_idx;%(y_u_idx*8/10 + y_s_idx*2/10); %y0_SE = 4.1937e-3; %35350 km from Earth
xdot0_SE = xdot_u_idx;
ydot0_SE = ydot0SE(1) %(ydot_u_idx*9/10 + ydot_s_idx*1/10);

S0_SE = [x0_SE; y0_SE; xdot0_SE; ydot0_SE];

figure
plot(y0_SE,ydot0_SE,'go')

% Insert into SE model (integrate Earth escape backwards in time from time patch in SE)
%Load
load('pos_SE','pos_SE')
load('states_po_SE','states_po_SE')
load('times_po_SE','times_po_SE')

% Integrate backwards
time_patch_SE = times_s_EM_SE{idx_EM}(end); %Time at patching point in SE (from Moon capture leg)
deltat = -4;%%Why -4
et0 = time_patch_SE;
prnt_out_dt = 0.01;

stop_fun = @(et,states_aux)earth_SE_crossing_detection_dec(et,states_aux,params_SE);

[SF_SE_SE_mid, etf_SE_SE_mid, states_SE_SE_mid, times_SE_SE_mid] = PCR3BP_propagator (S0_SE, params_SE,...
    et0, deltat, prnt_out_dt, stop_fun);

r = sqrt((SF_SE_SE_mid(1)-params_SE.mu1)^2 + SF_SE_SE_mid(2)^2)*L_SE
rdot = SF_SE_SE_mid(3).*cos(params_SE.phiA_tgt)-SF_SE_SE_mid(4).*sin(120*pi/180)*L_SE/T_SE

% Integrate forward 
deltat = (times_SE_SE_mid(1)-times_SE_SE_mid(end));
et0 = etf_SE_SE_mid;
prnt_out_dt = 0.001;

stop_fun = 'none';%@(et,states_aux)earth_SE_crossing_detection_inc(et,states_aux,params_SE);

[SF_SE_SE_result, etf_SE_SE_result, states_SE_SE_result, times_SE_SE_result] = PCR3BP_propagator (SF_SE_SE_mid, params_SE,...
    et0, deltat, prnt_out_dt, stop_fun);

%Save variables
save('states_SE_SE_result','states_SE_SE_result')
save('times_SE_SE_result','times_SE_SE_result')
save('SF_SE_SE_result','SF_SE_SE_result')
save('etf_SE_SE_result','etf_SE_SE_result')

%Plot
plot_Earth_scape_leg;

%% 3. Insert into EM model
%Load
load('pos_EM','pos_EM')
load('states_po_EM','states_po_EM')
load('times_po_EM','times_po_EM')

%Initial Conditions for integration => SO_EM_SE
 %ff = 3.2/8;%3/8;%2.7
 %xdot0_EM_SE = xdot_u_idx*(1-ff) + xdot_s_idx*ff;%(0.6*xdot_s_idx+0.4*xdot_u_idx); %0.5*(xdot_s_idx+xdot_s_idx2);
xdot0_EM_SE = 0.5*(xdot_s_idx+xdot_u_idx);
ydot0_EM_SE = ydot0_SE;
y0_EM_SE = y0_SE;
x0_EM_SE = x0_SE;

S0_EM_SE = [x0_EM_SE;y0_EM_SE;xdot0_EM_SE;ydot0_EM_SE];
S0_EM = SE2EM(S0_EM_SE, time_patch_SE, mu_EM, mu_SE, phi_EM_0, L_EM, L_SE, T_EM, T_SE);
S0_EM_SE_new = EM2SE(S0_EM, time_patch_SE*T_SE/T_EM, mu_EM, mu_SE, phi_EM_0, L_EM, L_SE, T_EM, T_SE);

%Integrate
et0 = time_patch_SE*T_SE/T_EM;
deltat = -time_patch_SE*T_SE/T_EM;%-time_patch_SE*T_SE/T_EM*1.2;
prnt_out_dt = 0.01;

stop_fun = @(et,states_aux)Moon_EM_crossing_detection(et,states_aux,params_EM);

[SF_EM_EM_result, etf_EM_EM_result, states_EM_EM_result, times_EM_EM_result] = PCR3BP_propagator (S0_EM, params_EM,...
    et0, deltat, prnt_out_dt, stop_fun);

xfEM = SF_EM_EM_result(1);
ydotfEM = SF_EM_EM_result(4);
yfEM  = SF_EM_EM_result(2);
xdotfEM = SF_EM_EM_result(3); 
rfEM = sqrt((xfEM-1+mu_EM)^2 + yfEM^2);
vfEM = sqrt((xdotfEM-yfEM)^2 + (ydotfEM + xfEM + mu_EM-1)^2);
veEM = sqrt(2*mu_EM/rfEM);
veEM-vfEM

%Save
save('states_EM_EM_result','states_EM_EM_result')
save('times_EM_EM_result','times_EM_EM_result')
save('SF_EM_EM_result','SF_EM_EM_result')
save('etf_EM_EM_result','etf_EM_EM_result')

% Transform Moon capture leg to SE system

for jj = 1:size(states_EM_EM_result,2)
    states_EM_SE_result(:,jj) = EM2SE(states_EM_EM_result(:,jj), times_EM_EM_result(jj), mu_EM, mu_SE, phi_EM_0, L_EM, L_SE, T_EM, T_SE);
    times_EM_SE_result(jj) = times_EM_EM_result(jj)*T_EM/T_SE;
end
SF_EM_SE_result = EM2SE(SF_EM_EM_result, etf_EM_EM_result, mu_EM, mu_SE, phi_EM_0, L_EM, L_SE, T_EM, T_SE);
etf_EM_SE_result = etf_EM_EM_result*T_SE/T_EM;
    
save('states_EM_SE_result','states_EM_SE_result')
save('times_EM_SE_result','times_EM_SE_result')
save('SF_EM_SE_result','SF_EM_SE_result')
save('etf_EM_SE_result','etf_EM_SE_result')

plot_Moon_capture_leg;

%% 4. Plot complete trajectory

states_SE_complete = [states_SE_SE_result states_EM_SE_result];
times_SE_complete = [times_SE_SE_result times_EM_SE_result];
states_EM_complete = zeros(size(states_SE_complete));
times_EM_complete = zeros(size(times_SE_complete));

for jj = 1:size(states_SE_complete,2)
    states_EM_complete(:,jj) = SE2EM(states_SE_complete(:,jj), times_SE_complete(jj), mu_EM, mu_SE, phi_EM_0, L_EM, L_SE, T_EM, T_SE);
    times_EM_complete(jj) = times_SE_complete(jj)*T_SE/T_EM;
end

save('states_EM_complete','states_EM_complete')
save('times_EM_complete','times_EM_complete')
save('states_SE_complete','states_SE_complete')
save('times_SE_complete','times_SE_complete')


plot_complete_trajectory;
    %Plot in SE
    %Plot in EM
    %Plot in inertial Earth centered
    %Plot in Inertial Moon centered
        