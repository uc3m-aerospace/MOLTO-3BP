function state_EM_SE = EM2SE(state_EM, times_EM, mu_EM, mu_SE, phi_EM_0, L_EM, L_SE, T_EM, T_SE)

%% Transform to inertial centered on the Earth (primary mass)
mu1_EM = 1-mu_EM;
mu2_EM = mu_EM;

state_EM_0 = [-mu2_EM; 0; 0; 0]; %Primary mass position

phi_tEM = times_EM+phi_EM_0;

c = cos(phi_tEM);
s = sin(phi_tEM);

R_11 = [c -s; s c];
R_22 = R_11;
R_21 = [-s -c; c -s];

R = [R_11 zeros(size(R_11)); R_21 R_22];

state_EM_in = R*(state_EM-state_EM_0);

%% Transform units

state_EM_in_dim(1:2,:) = state_EM_in(1:2,:)*L_EM/L_SE;
state_EM_in_dim(3:4,:) = state_EM_in(3:4,:)*L_EM/L_SE/(T_EM/T_SE);
times_EM_dim = times_EM*T_EM/T_SE;

%% Transform to SE rotating reference frame

c = cos(times_EM_dim);
s = sin(times_EM_dim);
R_11 = [c -s; s c];
R_22 = R_11;
R_21 = [-s -c; c -s];
R = [R_11 zeros(size(R_11)); R_21 R_22];

state_EM_SE = inv(R)*(state_EM_in_dim)+[(1-mu_SE);0;0;0];

end