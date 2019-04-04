%Plot Earth-scape leg
figure
plot(states_SE_SE_result(1,:),states_SE_SE_result(2,:),'b','linewidth',1.5)
hold on
plot(mu1_SE,0,'ko','markerfacecolor','k','markersize',4)%interp2(X,Y,U,-mu2,0),'ro')
text(mu1_SE-0.17,0,'$m_E$','interpreter','latex','fontsize',15)
plot(pos_SE(2,1),pos_SE(2,2),'bo','markerfacecolor','b','markersize',4)
text(pos_SE(2,1)+0.05,pos_SE(2,2),'$L_2$','interpreter','latex','fontsize',15)
title('Earth escape leg','interpreter','latex')
plot(states_po_SE(1,:),states_po_SE(2,:),'k','Linewidth',1.2)
grid on
xlabel('$X_{SE}$','interpreter','latex','fontsize',15)
ylabel('$Y_{SE}$','interpreter','latex','fontsize',15)
%set(gca,'ticklabelinterpreter','latex','fontsize',15)
%savefig('./Koon_v4/Initial_Guess_0304/Earth_escape_leg.fig')