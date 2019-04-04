function  plot_Manifolds(mu1,mu2,pos,states_po, states_s,SF_s,states_u,titstring)
%Inputs: load_variables, states_s_EM, states_u_EM, states_po_EM, pos_EM

figure
for ii = 1:size(states_s,2)
    AAA = size(states_s{ii}(1,:),2);
    SSline = plot(states_s{ii}(1,round(linspace(1,AAA,AAA/100))),states_s{ii}(2,round(linspace(1,AAA,AAA/100))),'b'); SSline.Color(4) = 0.3;
    hold on
    AAA = size(states_u{ii}(1,:),2);
    SUline = plot(states_u{ii}(1,round(linspace(1,AAA,AAA/3))),states_u{ii}(2,round(linspace(1,AAA,AAA/3))),'r'); SUline.Color(4) = 0.3;
end
plot(mu1,0,'ko','markerfacecolor','k','markersize',4)%interp2(X,Y,U,-mu2,0),'ro')
text(mu1+0.001,-0.004,'$m_M$','interpreter','latex','fontsize',18)
plot(pos(2,1),pos(2,2),'ko','markerfacecolor','k','markersize',4)
text(pos(2,1)-0.005,pos(2,2)-0.008,'$L_2$','interpreter','latex','fontsize',18)
plot(states_po(1,:),states_po(2,:),'k','Linewidth',2)
% plot(-mu2,0,'ko','markerfacecolor','k','markersize',4)%interp2(X,Y,U,-mu2,0),'bo')
% text(-mu2-0.17,0,'$m_E$','interpreter','latex','fontsize',15)
grid on
xlabel(strcat(['X-' titstring]),'interpreter','latex','fontsize',18)
ylabel(strcat(['Y-' titstring]),'interpreter','latex','fontsize',18)
axis([0.98 1.3 -0.11 0.11])
title('Manifolds','interpreter','latex')
%print('./Koon_v4/Data_Periodic_Orbits/Plots_Manifolds/EM_manifolds','-painters','-depsc2')
%set(gca,'ticklabelinterpreter','latex','fontsize',16)

figure
for ii = 1:size(states_u,2)
        plot(states_u{ii}(1,:),states_u{ii}(2,:),'r')
        hold on
        %%plot(states_s{ii}(1,:),states_s{ii}(2,:),'b')ARREGLARLOOO
        %%hold on
end
plot(-mu2,0,'bo','markerfacecolor','k','markersize',2)%interp2(X,Y,U,-mu2,0),'bo')
 text(-mu2,0,'$m_E$','interpreter','latex')
plot(mu1,0,'ro','markerfacecolor','k','markersize',2)%interp2(X,Y,U,-mu2,0),'ro')
 text(mu1,0,'$m_M$','interpreter','latex')
plot(pos(2,1),pos(2,2),'ko','markerfacecolor','k','markersize',2)
 text(pos(2,1),pos(2,2),'$L_2$','interpreter','latex')
title(strcat(['Stable and Unstable manifolds' titstring]),'interpreter','latex')
plot(states_po(1,:),states_po(2,:),'k','Linewidth',1.2)
%savefig('EM_manifolds.fig')
    
figure
plot(-mu2,0,'bo','markerfacecolor','k','markersize',1)%interp2(X,Y,U,-mu2,0),'bo')
text(-mu2,0,'$m_E$','interpreter','latex')
plot(mu1,0,'ro','markerfacecolor','k','markersize',1)%interp2(X,Y,U,-mu2,0),'ro')
text(mu1,0,'$m_M$','interpreter','latex')
plot(pos(1:2,1),pos(1:2,2),'ko','markerfacecolor','k','markersize',1)
title(strcat(['Stable manifolds' titstring]),'interpreter','latex')
       
figure
plot(SF_s(2,:),SF_s(4,:),'o')
hold on
xlabel('$\bar{y}$','interpreter','latex')
ylabel('$\dot{y}$','interpreter','latex')
title('Phase space in Poincare section','interpreter','latex')
%savefig('EM_phspc.fig')

end