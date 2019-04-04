%Plot Moon Capture leg
figure
plot(states_EM_EM_result(1,:),states_EM_EM_result(2,:),'b','linewidth',1.5)
hold on
plot(-mu2_EM,0,'ko','markerfacecolor','k','markersize',4)%interp2(X,Y,U,-mu2,0),'bo')
text(-mu2_EM-0.17,0,'$m_E$','interpreter','latex','fontsize',15)
plot(mu1_EM,0,'ko','markerfacecolor','k','markersize',4)%interp2(X,Y,U,-mu2,0),'ro')
text(mu1_EM-0.17,0,'$m_M$','interpreter','latex','fontsize',15)
plot(pos_EM(2,1),pos_EM(2,2),'bo','markerfacecolor','b','markersize',4)
text(pos_EM(2,1)+0.05,pos_EM(2,2),'$L_2$','interpreter','latex','fontsize',15)
title('Moon capture leg','interpreter','latex')
plot(states_po_EM(1,:),states_po_EM(2,:),'k','Linewidth',1.2)
grid on
xlabel('$X_{EM}$','interpreter','latex','fontsize',15)
ylabel('$Y_{EM}$','interpreter','latex','fontsize',15)
%set(gca,'ticklabelinterpreter','latex','fontsize',15)
%savefig('./Koon_v4/Initial_Guess_0304/Moon_capture_leg_EM.fig')