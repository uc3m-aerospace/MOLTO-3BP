% phi_A = 90*pi/180;
% 
% nn_SE_s_min = 1;
% exit_phiA = 0;
% nn_SE_s_inf = 1;
% 
% for ii = 1:size(states_s_SE,2)
%     if states_s_SE{ii}(1,end) < pos_SE(2,1)
%         for jj = 1:size(states_s_SE{ii},2)
%             
%             % Stable manifold EM system (yplus part)
%             phi_A_t = atan2(-states_s_SE{ii}(2,jj),(states_s_SE{ii}(1,jj)-params_SE.mu1));
%             
%             if (phi_A_t < phi_A) || states_s_SE{ii}(2,jj) > 0
%                 states_s_SE_phiA{nn_SE_s_min}(:,jj) = states_s_SE{ii}(:,jj);
%                 times_s_SE_phiA{nn_SE_s_min}(jj) = times_s_SE{ii}(jj);
%                 exit_phiA = 1;
%                 %                 elseif phi_A_t > phi_A && exit_phiA == 1
%                 %                     SF_s_SE_phiA(:,nn_SE_s_min) = states_s_SE_phiA{nn_SE_s_min}(:,end);
%             end
%         end
%         SF_s_SE_phiA(:,nn_SE_s_min) = states_s_SE_phiA{nn_SE_s_min}(:,end);
%         
%         nn_SE_s_min = nn_SE_s_min+1;
%         
%    elseif states_s_SE{ii}(1,end) > pos_SE(2,1)
%         states_s_SE_inf{nn_SE_s_inf} = states_s_SE{ii};
%         nn_SE_s_inf = nn_SE_s_inf+1;
%     end
% end
% % save(strcat(['./Koon_v4/Data_Periodic_Orbits/Plots_Manifolds/states_s_SE_phiA_' num2str(phi_A*180/pi)]),'states_s_SE_phiA')
% % save(strcat(['./Koon_v4/Data_Periodic_Orbits/Plots_Manifolds/times_s_SE_phiA_' num2str(phi_A*180/pi)]),'times_s_SE_phiA')
% % save(strcat(['./Koon_v4/Data_Periodic_Orbits/Plots_Manifolds/SF_s_SE_phiA_' num2str(phi_A*180/pi)]),'SF_s_SE_phiA')
% % save(strcat(['./Koon_v4/Data_Periodic_Orbits/Plots_Manifolds/states_s_SE_inf_phiA_']),'states_s_SE_inf')
% 
% phi_B = 90*pi/180;
% 
% nn_SE_u_plus = 1;
% nn_SE_u_inf = 1;
% exit_phiB = 0;
% for ii = 1:size(states_u_SE,2)
%     if states_u_SE{ii}(1,end) < pos_SE(2,1)
%         for jj = 1:size(states_u_SE{ii},2)
%             
%             % Stable manifold EM system (yplus part)
%             phi_B_t = atan2(states_u_SE{ii}(2,jj),(states_u_SE{ii}(1,jj)-params_SE.mu1));
%             
%             if (phi_B_t < phi_B) || states_u_SE{ii}(2,jj) < 0
%                 states_u_SE_phiB{nn_SE_u_plus}(:,jj) = states_u_SE{ii}(:,jj);
%                 times_u_SE_phiB{nn_SE_u_plus}(jj) = times_u_SE{ii}(jj);
%                 exit_phiB = 1;
%             end
%         end
%         SF_u_SE_phiB(:,nn_SE_u_plus) = states_u_SE_phiB{nn_SE_u_plus}(:,end);
%         
%         nn_SE_u_plus = nn_SE_u_plus+1;
%         
%     else
%         states_u_SE_inf{nn_SE_u_inf} = states_u_SE{ii};
%         nn_SE_u_inf = nn_SE_u_inf+1;
%     end
% end
% % save(strcat(['./Koon_v4/Data_Periodic_Orbits/Plots_Manifolds/states_u_SE_phiB_' num2str(phi_B*180/pi)]),'states_u_SE_phiB')
% % save(strcat(['./Koon_v4/Data_Periodic_Orbits/Plots_Manifolds/times_u_SE_phiB_' num2str(phi_B*180/pi)]),'times_u_SE_phiB')
% % save(strcat(['./Koon_v4/Data_Periodic_Orbits/Plots_Manifolds/SF_u_SE_phiB_' num2str(phi_B*180/pi)]),'SF_u_SE_phiB')


figure
for ii = 1:nn_SE_s_min-1
    %figure(3)
    AAA = size(states_s_SE_phiA{ii}(1,:),2);
    SAline = plot(states_s_SE_phiA{ii}(1,round(linspace(1,AAA,AAA/300))),states_s_SE_phiA{ii}(2,round(linspace(1,AAA,AAA/300))),'b'); SAline.Color(4) = 0.3;
    hold on
    AAA = size(states_s_SE_inf{ii}(1,:),2);
    Sinfline = plot(states_s_SE_inf{ii}(1,round(linspace(1,AAA,AAA/300))),states_s_SE_inf{ii}(2,round(linspace(1,AAA,AAA/300))),'b'); Sinfline.Color(4) = 0.3;
end

for ii = 1:nn_SE_u_plus-1
    %figure(3)
    AAA = size(states_u_SE_phiB{ii}(1,:),2);
    SAline = plot(states_u_SE_phiB{ii}(1,round(linspace(1,AAA,AAA/30))),states_u_SE_phiB{ii}(2,round(linspace(1,AAA,AAA/30))),'r'); SAline.Color(4) = 0.3;
    hold on
    AAA = size(states_u_SE_inf{ii}(1,:),2);
    Sinfline = plot(states_u_SE_inf{ii}(1,round(linspace(1,AAA,AAA/30))),states_u_SE_inf{ii}(2,round(linspace(1,AAA,AAA/30))),'r'); Sinfline.Color(4) = 0.3;
end

figure
plot(states_po_SE(1,:),states_po_SE(2,:),'k','Linewidth',2)
grid on
plot(mu1_SE,0,'ko','markerfacecolor','k','markersize',4)%interp2(X,Y,U,-mu2,0),'ro')
text(mu1_SE-0.0005,-0.0005,'$m_E$','interpreter','latex','fontsize',18,'fontweight','bold')
plot(pos_SE(2,1),pos_SE(2,2),'ko','markerfacecolor','k','markersize',4)
text(pos_SE(2,1)-0.0005,pos_SE(2,2)-0.0005,'$L_2$','interpreter','latex','fontsize',18,'fontweight','bold')
xlabel('$X_{SE}$','interpreter','latex')
ylabel('$Y_{SE}$','interpreter','latex')
%set(gca,'ticklabelinterpreter','latex','fontsize',16)
axis([0.999 1.02 -0.008 0.008])
%print('./Koon_v4/Data_Periodic_Orbits/Plots_Manifolds/SE_manifolds','-painters','-depsc2')
% savefig('./Koon_v4/Data_Periodic_Orbits/Plots_Manifolds/SE_manifolds.fig')