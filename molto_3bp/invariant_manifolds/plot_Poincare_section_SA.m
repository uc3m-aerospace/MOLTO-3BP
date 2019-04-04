function plot_Poincare_section_SA(state, fitfactor, phiA, titstring,xE,R)


x = state(1,:)-xE;
y = state(2,:);
r = sqrt(x.^2+y.^2);
xdot = state(3,:);
ydot = state(4,:);
rdot = xdot.*cos(phiA)-ydot.*sin(phiA);

figure
k = boundary(y',ydot',fitfactor);
plot(y,ydot,'o','color','b','markersize',3,'markerfacecolor','b')
hold on
xlabel('$Y_{SE}$','interpreter','latex','fontsize',15)
ylabel('$\dot{Y}_{SE}$','interpreter','latex','fontsize',15)
title(titstring,'interpreter','latex','fontsize',15)
%set(gca,'ticklabelinterpreter','latex','fontsize',15)
ylim=get(gca,'ylim');
xlim=get(gca,'xlim');
grid on
plot([-R -R],[ylim(1) ylim(2)],'k--')
%savefig(strcat(['poincare_SA_phiA_' num2str(phiA*180/pi) '.fig']))
%print(strcat(['poincare_SA_phiA_' num2str(phiA*180/pi)]),'-depsc2')
end