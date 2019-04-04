%POST-PROCESSING OF MANIFOLDS
%The trajectory will depart from GEO at SA to the SE stable manifold, then twist
%towards the SE unstable manifold (Earth scape leg), and at SB, impulse towards EM stable
%manifold (Moon capture leg)

%It is divided in the following sections:

    %%% 0. Initialize General variables
    %%% 1. MOON-LEG: Transform Poincaré sections to Bicircular Reference frame
    %%% 2. GET SA POINCARÉ SECTION (SE STABLE MANIFOLD)
    %%% 3. GET SB POINCARÉ SECTION (UNSTABLE MANIFOLD)
    %%% 4. 4. GET POINCARÉ SECTION PATH
        

%%0. Initialize General variables
close all
clear all
clc 
%Mice package, Gravitational constants, Earth-Moon, and Sun-Earth constants
load_variables;

%Initial angle formed by x-axis of SE and EM system
phi_EM_0 = 225*pi/180;
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

%%  1. Transform Poincaré sections to Bicircular Reference frame

    %1.1 Load stable manifolds of Earh-Moon system
    %load(filename,variables) loads the specified variables from the MAT-file, filename.
    %load('./"FOlder route"/SF_s_EM','SF_s_EM')
    
    load(strcat(['states_s_EM_phi0_' num2str(phi_EM_0*180/pi)]),'states_s_EM')
    load(strcat(['times_s_EM_phi0_'  num2str(phi_EM_0*180/pi)]),'times_s_EM')
    load(strcat(['SF_s_EM_phi0_'  num2str(phi_EM_0*180/pi)]),'SF_s_EM')

    %% 1.2 This script computes the Moon leg
%     bcrf = 1;
%     if bcrf == 1 
%         load(strcat(['states_s_EM_SE' num2str(phi_EM_0*180/pi)]),'states_s_EM_SE')
%         load(strcat(['times_s_EM_SE' num2str(phi_EM_0*180/pi)]),'times_s_EM_SE')
%         load(strcat(['SF_s_EM_SE' num2str(phi_EM_0*180/pi)]),'SF_s_EM_SE')
%     else
        Poincare_to_BCRF;
        %1.3 Outputs:
        save(strcat(['states_s_EM_SE' num2str(phi_EM_0*180/pi)]),'states_s_EM_SE')
        save(strcat(['times_s_EM_SE' num2str(phi_EM_0*180/pi)]),'times_s_EM_SE')
        save(strcat(['SF_s_EM_SE' num2str(phi_EM_0*180/pi)]),'SF_s_EM_SE')
%     end
    

    %% 1.3 Plots
        %a .Exterior stable Earth-Moon manifold
        %b. Phase space at poincare section SB
    fitfactor = 0.6;
    titstring = strcat(['Phase space at Poincare section $S_B$, $\varphi_0$ = ' num2str(phi_EM_0*180/pi) '$^\circ$']); 
    
    plot_Poincare_section_S0%(SF_s_EM_SE, fitfactor, phi_EM_0, titstring,params_SE.mu1, R_SE)
  
      
%%  2. GET SA POINCARÉ SECTION (SE STABLE MANIFOLD)

    %2.1 Load stable manifolds of Sun-Earh system, the position of L points
    %and states of periodic orbits
    load(strcat(['states_s_SE_phi0_' num2str(phi_EM_0*180/pi)]),'states_s_SE')
    load(strcat(['times_s_SE_phi0_'  num2str(phi_EM_0*180/pi)]),'times_s_SE')
    load(strcat(['SF_s_SE_phi0_'  num2str(phi_EM_0*180/pi)]),'SF_s_SE')
    load('pos_SE','pos_SE')
    load('states_po_SE','states_po_SE')
    
    %% 2.2 This script computes phiA
    phi_A = 90*pi/180;%Why 150
    
%     Sa=0;
%     if Sa==1
%          load(strcat(['states_s_SE_phiA_' num2str(phi_A*180/pi)]),'states_s_SE_phiA')
%          load(strcat(['times_s_SE_phiA_' num2str(phi_A*180/pi)]),'times_s_SE_phiA')
%          load(strcat(['SF_s_SE_phiA_' num2str(phi_A*180/pi)]),'SF_s_SE_phiA')
%     else
    SA_poincare_section;
    
    %2.3 Outputs
    save(strcat(['states_s_SE_phiA_' num2str(phi_A*180/pi)]),'states_s_SE_phiA')
    save(strcat(['times_s_SE_phiA_' num2str(phi_A*180/pi)]),'times_s_SE_phiA')
    save(strcat(['SF_s_SE_phiA_' num2str(phi_A*180/pi)]),'SF_s_SE_phiA')
%     end 
    
   % Get Poincaré section on S_A with varying phi_A (SE stable manifold)

   
   %% 2.4 plots
   fitfactor = 0.6;
   titstring = strcat(['Poincare section $S_A$, $\varphi_A$ = ' num2str(phi_A*180/pi) '$^\circ$']);
        
   plot_Poincare_section_SA(SF_s_SE_phiA, fitfactor, phi_A, titstring,params_SE.mu1, R_SE)
   
%% 3. GET SB POINCARÉ SECTION (UNSTABLE MANIFOLD) 

    %3.1 Load unstable manifolds of Sun-Earh system 
    load(strcat(['states_u_SE_phi0_' num2str(phi_EM_0*180/pi)]),'states_u_SE')
    load(strcat(['times_u_SE_phi0_'  num2str(phi_EM_0*180/pi)]),'times_u_SE')
    load(strcat(['SF_u_SE_phi0_'  num2str(phi_EM_0*180/pi)]),'SF_u_SE')
    load('pos_SE','pos_SE')
    load('states_po_SE','states_po_SE')
    load(strcat(['states_s_SE_phi0_' num2str(phi_EM_0*180/pi)]),'states_s_SE')
        
    %% 3.2 This script computes phiB
    phi_B = 90*pi/180;
%     Sb=1;
%     if Sb==0
%         load(strcat(['states_u_SE_phiB_' num2str(phi_B*180/pi)]),'states_u_SE_phiB')
%         load(strcat(['times_u_SE_phiB_' num2str(phi_B*180/pi)]),'times_u_SE_phiB')
%         load(strcat(['SF_u_SE_phiB_' num2str(phi_B*180/pi)]),'SF_u_SE_phiB')
%     else
        SB_poincare_section;
        %3.3 Outputs
        save(strcat(['states_u_SE_phiB_' num2str(phi_B*180/pi)]),'states_u_SE_phiB')
        save(strcat(['times_u_SE_phiB_' num2str(phi_B*180/pi)]),'times_u_SE_phiB')
        save(strcat(['SF_u_SE_phiB_' num2str(phi_B*180/pi)]),'SF_u_SE_phiB')
%     end
    %% 3.4 Plots
    fitfactor = 1;
    titstring = strcat(['Poincare section $S_B$, $\varphi_B$ = ' num2str(phi_B*180/pi) '$^\circ$']);
    
    plot_Poincare_section_SB(SF_u_SE_phiB, fitfactor, phi_B, titstring,params_SE.mu1, R_SE)
%%  4. GET POINCARÉ SECTION PATH(...)
%Inputs
%load(strcat(['states_s_EM_phi0_' num2str(phi_EM_0*180/pi)]),'states_s_EM_SE')
%load(strcat(['times_s_EM_phi0_'  num2str(phi_EM_0*180/pi)]),'times_s_EM_SE')
%load(strcat(['SF_s_EM_phi0_'  num2str(phi_EM_0*180/pi)]),'SF_u_SE_phiB')

titstring = strcat(['Phase space at Poincare section $S_B$, $\varphi_0$ = ' num2str(phi_EM_0*180/pi) '$^\circ$, $\varphi_B$ = ' num2str(90) '$^\circ$']);
fitfactor = 0;

plot_Poincare_section_patch(SF_u_SE_phiB, SF_s_EM_SE, fitfactor, phi_B, titstring, mu1_SE, R_SE)
   
%% 5. 
%plot_raro;
   
   