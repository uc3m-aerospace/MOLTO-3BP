function plot_Poincare_section_SA(state_SE, state_EM, fitfactor, phiB, titstring,xE,R)

x_SE = state_SE(1,:)-xE;
y_SE = state_SE(2,:);
xdot_SE = state_SE(3,:);
ydot_SE = state_SE(4,:);

figure
k = boundary(y_SE',ydot_SE',fitfactor);
plot(y_SE,ydot_SE,'o','color','r','markersize',2,'markerfacecolor','r')
hold on
% plot(y_SE(k),ydot_SE(k),'color','r')
%h = fill(y_SE,ydot_SE,'r','edgecolor','none');
%set(h,'facealpha',.5)

xlabel('$\bar{y}$','interpreter','latex','fontsize',15)
ylabel('$\dot{y}$','interpreter','latex','fontsize',15)
title(titstring,'interpreter','latex','fontsize',15)
%set(gca,'ticklabelinterpreter','latex','fontsize',15)
ylim=get(gca,'ylim');
xlim=get(gca,'xlim');
grid on
plot([R R],[ylim(1) ylim(2)],'k--')
%savefig(strcat(['poincare_SB_phiB_' num2str(phiB*180/pi) '.fig']))

x_EM = state_EM(1,:)-xE;
y_EM = state_EM(2,:);
ydot_EM = state_EM(4,:);

k = boundary(y_EM',ydot_EM',fitfactor);
plot(y_EM,ydot_EM,'o','color','b','markersize',2,'markerfacecolor','b')
hold on
%plot(y_EM(k),ydot_EM(k),'color','b')
% h = fill(y_EM(k),ydot_EM(k),'b','edgecolor','none');
% set(h,'facealpha',.5)

end