%function plot_Poincare_section_S0(state, fitfactor, phi_EM_0, titstring,xE,R)

% %Plots    
figure    
for ii = 1:size(states_s_EM,2)-1%nn_EM_s_plus-1
        %figure
        SCline = plot(states_s_EM_SE{ii}(1,:),states_s_EM_SE{ii}(2,:),'b');  SCline.Color(4) = 0.3;
        hold on
end
    
states_Moon_0 = EM2SE([params_EM.mu1;0;0;0],0, params_EM.mu_EM, params_EM.mu_SE, params_EM.phi_EM_0, params_EM.L_EM, params_EM.L_SE, params_EM.T_EM, params_EM.T_SE);
moon_orbit_SE = circle(mu1_SE,0,L_EM/L_SE);
    
%figure
grid on
axis equal
plot(mu1_SE,0,'ko','markerfacecolor','k','markersize',6)%interp2(X,Y,U,-mu2,0),'ro')
text(mu1_SE*1.0001,0,'$m_E$','interpreter','latex','fontsize',15,'fontweight','bold')
plot(states_Moon_0(1),states_Moon_0(2),'ko','markerfacecolor','k','markersize',6)
text(states_Moon_0(1)*0.9993,states_Moon_0(2),'$m_{M,0}$','interpreter','latex','fontsize',15,'fontweight','bold')
title(strcat(['Earth-Moon stable manifold, $\varphi_0$ = ' num2str(phi_EM_0*180/pi) '$^\circ$']),'interpreter','latex','fontsize',15)
plot(moon_orbit_SE(:,1),moon_orbit_SE(:,2),'k--','linewidth',1)
xlabel('$\bar{x}$','interpreter','latex')
ylabel('$\bar{y}$','interpreter','latex')
ylim = get(gca,'ylim');
plot([mu1_SE,mu1_SE],[0 ylim(2)],'k','Linewidth',2)
text(mu1_SE*1.0001,ylim(2)*0.5,'$S_0$','interpreter','latex','fontsize',15,'fontweight','bold')
%set(gca,'ticklabelinterpreter','latex','fontsize',15)
% savefig(strcat(['EM_transfer_manifolds_phi0_' num2str(phi_EM_0*180/pi) '.fig']))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%state = SF_s_EM_SE;
%xE = params_SE.mu1;

x = SF_s_EM_SE(1,:)-params_SE.mu1;
y = SF_s_EM_SE(2,:);
ydot = SF_s_EM_SE(4,:);
%openfig('./Koon_v4/Data_Periodic_Orbits/SB_manifolds/poincare_SB_phiB_90.fig','new');
figure
k = boundary(y',ydot',fitfactor);
plot(y,ydot,'o','color','b','markersize',2,'markerfacecolor','b')
hold on
%plot(y(k),ydot(k),'color','b')
%h = fill(y(k),ydot(k),'b','edgecolor','none');
%set(h,'facealpha',.5)
xlabel('$Y_{SE}$','interpreter','latex','fontsize',18)
ylabel('$\dot{Y}_{SE}$','interpreter','latex','fontsize',18)
title(titstring,'interpreter','latex','fontsize',15)
%set(gca,'ticklabelinterpreter','latex','fontsize',16)
grid on
axis([0 5e-3 -0.12 0])
ylim=get(gca,'ylim');
xlim=get(gca,'xlim');
plot([R_SE R_SE],[ylim(1) ylim(2)],'k--')
%savefig(strcat(['./Koon_v4/Data_Periodic_Orbits/EM_manifolds/poincare_S0_phi0_' num2str(phi_EM_0*180/pi) '.fig']))
%print(strcat(['./Koon_v4/Data_Periodic_Orbits/EM_manifolds/poincare_S0_phi0_' num2str(phi_EM_0*180/pi)]),'-painters','-depsc2')
%end