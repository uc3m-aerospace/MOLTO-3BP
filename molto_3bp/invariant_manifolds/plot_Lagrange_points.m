function  plot_Lagrange_points(mu1,mu2, pos,titstring)

figure
plot(-mu2,0,'ko','markerfacecolor','k','markersize',5)
hold on
text(-mu2-0.05,-0.07,'$m_E$','interpreter','latex','fontsize',18)
plot(mu1,0,'ko','markerfacecolor','k','markersize',5)
text(mu1-0.05,-0.07,'$m_M$','interpreter','latex','fontsize',18)
plot(pos(1,1),pos(1,2),'bo','markerfacecolor','b','markersize',5)
text(pos(1,1)-0.15,pos(1,2),'$L_1$','interpreter','latex','fontsize',18)
plot(pos(2,1),pos(2,2),'bo','markerfacecolor','b','markersize',5)
text(pos(2,1)+0.05,pos(2,2),'$L_2$','interpreter','latex','fontsize',18)
plot(pos(3,1),pos(3,2),'bo','markerfacecolor','b','markersize',5)
text(pos(3,1)+0.05,pos(3,2),'$L_3$','interpreter','latex','fontsize',18)
plot(pos(4,1),pos(4,2),'ro','markerfacecolor','r','markersize',5)
text(pos(4,1)+0.05,pos(4,2),'$L_4$','interpreter','latex','fontsize',18)
plot(pos(5,1),pos(5,2),'ro','markerfacecolor','r','markersize',5)
text(pos(5,1)+0.05,pos(5,2),'$L_5$','interpreter','latex','fontsize',18)
plot([-mu2 pos(4,1)],[0 pos(4,2)],'--','color',[0.5 0.5 0.5])
plot([pos(4,1) mu1],[pos(4,2) 0],'--','color',[0.5 0.5 0.5])
plot([-mu2 pos(5,1)],[0 pos(5,2)],'--','color',[0.5 0.5 0.5])
plot([pos(5,1) mu1],[pos(5,2) 0],'--','color',[0.5 0.5 0.5])

moon_orbit = circle(-mu2,0,1);
plot(moon_orbit(:,1),moon_orbit(:,2),'--','color',[0.5 0.5 0.5])
axis equal
title(strcat(['Lagrange points for the ' titstring ' system']),'interpreter','latex','fontsize',15)
xlabel(strcat(['X-' titstring]),'interpreter','latex','fontsize',15)
ylabel(strcat(['Y-' titstring]),'interpreter','latex','fontsize',15)
grid on
axis([-1.4 1.4 -1.2 1.2])
%set(gca,'ticklabelinterpreter','latex','fontsize',15)
%print('lagrange_EM','-depsc2')
end