function plot_Poincare_section_SA(state, fitfactor, phiB, titstring,xE,R)


x = state(1,:)-xE;
y = state(2,:);
xdot = state(3,:);
ydot = state(4,:);

figure
k = boundary(y',ydot',fitfactor);
plot(y,ydot,'o','color','r','markersize',2,'markerfacecolor','b')
hold on
%plot(y(k),ydot(k),'color','r')
%h = fill(y(k),ydot(k),'r','edgecolor','none');
%set(h,'facealpha',.5)
xlabel('$\bar{y}$','interpreter','latex','fontsize',15)
ylabel('$\dot{y}$','interpreter','latex','fontsize',15)
title(titstring,'interpreter','latex','fontsize',15)
set(gca,'ticklabelinterpreter','latex','fontsize',15)
ylim=get(gca,'ylim');
xlim=get(gca,'xlim');
grid on
plot([R R],[ylim(1) ylim(2)],'k--')
%savefig(strcat(['poincare_SB_phiB_' num2str(phiB*180/pi) '.fig']))
end