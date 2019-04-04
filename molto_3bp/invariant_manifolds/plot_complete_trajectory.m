%%%Plot complete trajectory
%% Plot in SE
states_Moon_f = EM2SE([params_EM.mu1;0;0;0],etf_EM_EM_result, params_EM.mu_EM, params_EM.mu_SE, params_EM.phi_EM_0, params_EM.L_EM, params_EM.L_SE, params_EM.T_EM, params_EM.T_SE);
states_Moon_patch = EM2SE([params_EM.mu1;0;0;0],time_patch_SE*T_SE/T_EM, params_EM.mu_EM, params_EM.mu_SE, params_EM.phi_EM_0, params_EM.L_EM, params_EM.L_SE, params_EM.T_EM, params_EM.T_SE);

moon_orbit_SE = circle(mu1_SE,0,L_EM/L_SE);

figure
plot(states_SE_complete(1,:),states_SE_complete(2,:),'k','Linewidth',2)
hold on
plot(mu1_SE,0,'ko','markerfacecolor','k','markersize',4)%interp2(X,Y,U,-mu2,0),'ro')
 text(mu1_SE-0.0009,0,'$m_E$','interpreter','latex','fontsize',15)
plot(pos_SE(2,1),pos_SE(2,2),'ko','markerfacecolor','k','markersize',4)
 text(pos_SE(2,1)+0.0003,pos_SE(2,2),'$L_2$','interpreter','latex','fontsize',15)
plot(states_po_SE(1,:),states_po_SE(2,:),'k','Linewidth',1.2)
grid on
plot(moon_orbit_SE(:,1),moon_orbit_SE(:,2),'k--','linewidth',1)
plot(states_Moon_f(1),states_Moon_f(2),'ko','markerfacecolor','k','markersize',6)
 text(states_Moon_f(1)-0.0005,states_Moon_f(2)-0.0005,'$m_{M,F}$','interpreter','latex','fontsize',15,'fontweight','bold')
plot(states_Moon_patch(1),states_Moon_patch(2),'ko','markerfacecolor','k','markersize',6)
 text(states_Moon_patch(1)-0.0005,states_Moon_patch(2)-0.0005,'$m_{M,P}$','interpreter','latex','fontsize',15,'fontweight','bold')
title('Transfer trajectory','interpreter','latex')
xlabel('$X_{SE}$','interpreter','latex','fontsize',15)
ylabel('$Y_{SE}$','interpreter','latex','fontsize',15)
%set(gca,'ticklabelinterpreter','latex','fontsize',15)
%savefig('./Koon_v4/Initial_Guess_GEO_0804/Complete_transfer_SE.fig')


load('states_s_SE_phiA_90','states_s_SE_phiA')
load('states_u_SE_phiB_90','states_u_SE_phiB')

figure
for ii = 1:size(states_s_SE_phiA,2)
    SAline = plot(states_s_SE_phiA{ii}(1,:),states_s_SE_phiA{ii}(2,:),'b'); SAline.Color(4) = 0.3;
end
for ii = 1:size(states_u_SE_phiB,2)
    SBline = plot(states_u_SE_phiB{ii}(1,:),states_u_SE_phiB{ii}(2,:),'r'); SBline.Color(4) = 0.3;
end

%% Plot in EM
figure
plot(states_EM_complete(1,:),states_EM_complete(2,:),'k','Linewidth',2)
hold on
plot(mu1_EM,0,'ko','markerfacecolor','k','markersize',4)%interp2(X,Y,U,-mu2,0),'ro')
 text(mu1_EM-0.0009,0,'$m_M$','interpreter','latex','fontsize',15)
plot(-mu2_EM,0,'ko','markerfacecolor','k','markersize',4)%interp2(X,Y,U,-mu2,0),'ro')
 text(-mu2_EM-0.0009,0,'$m_E$','interpreter','latex','fontsize',15)
plot(pos_EM(2,1),pos_EM(2,2),'ko','markerfacecolor','k','markersize',4)
 text(pos_EM(2,1)+0.0003,pos_EM(2,2),'$L_2$','interpreter','latex','fontsize',15)
title('Transfer trajectory','interpreter','latex')
grid on
xlabel('$X_{EM}$','interpreter','latex','fontsize',15)
ylabel('$Y_{EM}$','interpreter','latex','fontsize',15)
%set(gca,'ticklabelinterpreter','latex','fontsize',15)
%savefig('./Koon_v4/Initial_Guess_GEO_0804/Complete_transfer_EM.fig')

%% Plot in Inertial Earth centered

states_INE_complete = zeros(size(states_SE_complete));
times_INE_complete = zeros(size(times_SE_complete));

for jj = 1:size(states_SE_complete,2)
    states_INE_complete(:,jj) = EM2INE(states_EM_complete(:,jj), times_EM_complete(jj), mu_EM, mu_SE, phi_EM_0, L_EM, L_SE, T_EM, T_SE);
    times_INE_complete(jj) = times_EM_complete(jj)*T_SE/T_EM;
end

save('states_INE_complete','states_INE_complete')
save('times_INE_complete','times_INE_complete')

figure
plot(states_INE_complete(1,:),states_INE_complete(2,:),'k','Linewidth',2)
hold on
% plot(mu1_EM,0,'ko','markerfacecolor','k','markersize',4)%interp2(X,Y,U,-mu2,0),'ro')
% text(mu1_EM-0.0009,0,'$m_M$','interpreter','latex','fontsize',15)
plot(0,0,'ko','markerfacecolor','k','markersize',4)%interp2(X,Y,U,-mu2,0),'ro')
 text(0-0.0009,0,'$m_E$','interpreter','latex','fontsize',15)
plot(pos_EM(2,1),pos_EM(2,2),'ko','markerfacecolor','k','markersize',4)
 text(pos_EM(2,1)+0.0003,pos_EM(2,2),'$L_2$','interpreter','latex','fontsize',15)
title('Transfer trajectory','interpreter','latex')
grid on
xlabel('$X_{EM}$','interpreter','latex','fontsize',15)
ylabel('$Y_{EM}$','interpreter','latex','fontsize',15)
set(gca,'ticklabelinterpreter','latex','fontsize',15)
savefig('./Koon_v4/Initial_Guess_GEO_0804/Complete_transfer_INE.fig')

%% Plot in Inertial Moon centered

states_INM_complete = zeros(size(states_SE_complete));
times_INM_complete = zeros(size(times_SE_complete));

for jj = 1:size(states_SE_complete,2)
    states_INM_complete(:,jj) = EM2INM(states_EM_complete(:,jj), times_EM_complete(jj), mu_EM, mu_SE, phi_EM_0, L_EM, L_SE, T_EM, T_SE);
    times_INM_complete(jj) = times_EM_complete(jj)*T_SE/T_EM;
end

save('states_INM_complete','states_INM_complete')
save('times_INM_complete','times_INM_complete')

figure
plot(states_INM_complete(1,:),states_INM_complete(2,:),'k','Linewidth',2)
hold on
plot(0,0,'ko','markerfacecolor','k','markersize',4)%interp2(X,Y,U,-mu2,0),'ro')
text(0-0.0009,0,'$m_M$','interpreter','latex','fontsize',15)
% plot(-mu2_EM,0,'ko','markerfacecolor','k','markersize',4)%interp2(X,Y,U,-mu2,0),'ro')
% text(-mu2_EM-0.0009,0,'$m_E$','interpreter','latex','fontsize',15)
plot(pos_EM(2,1),pos_EM(2,2),'ko','markerfacecolor','k','markersize',4)
text(pos_EM(2,1)+0.0003,pos_EM(2,2),'$L_2$','interpreter','latex','fontsize',15)
title('Transfer trajectory','interpreter','latex')
grid on
xlabel('$X_{EM}$','interpreter','latex','fontsize',15)
ylabel('$Y_{EM}$','interpreter','latex','fontsize',15)
%set(gca,'ticklabelinterpreter','latex','fontsize',15)
%savefig('./Koon_v4/Initial_Guess_GEO_0804/Complete_transfer_INM.fig')


