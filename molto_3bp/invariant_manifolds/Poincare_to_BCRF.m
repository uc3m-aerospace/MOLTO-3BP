%Poincare to Bicircular Reference frame

% nn_EM_s_plus = 1;
%    for ii = 1:size(states_s_EM,2)
%         %SF_check = EM2SE(SF_s_EM(:,ii),times_s_EM{ii}(end),mu_EM,mu_SE,phi_EM_0,L_EM,L_SE,T_EM,T_SE);
%        % if (abs(SF_check(1) - mu1_SE) < 10^-4) && (SF_check(2) > 3e-3)
%             for jj = 1:size(states_s_EM{ii},2)
%                 states_s_EM_SE{nn_EM_s_plus}(:,jj) = EM2SE(states_s_EM{ii}(:,jj),times_s_EM{ii}(jj),mu_EM,mu_SE,phi_EM_0,L_EM,L_SE,T_EM,T_SE);
%                 times_s_EM_SE{nn_EM_s_plus}(:,jj) = times_s_EM{ii}(jj)*T_EM/T_SE;
%             end
%             SF_s_EM_SE(:,nn_EM_s_plus) = states_s_EM_SE{nn_EM_s_plus}(:,end);
%             nn_EM_s_plus = nn_EM_s_plus+1;
%        % end
%         ii
%     end

%nn_EM_s_plus = 1;
for ii = 1:size(states_s_EM,2)%size(states_s_EM,2) returns the number of column cells 
    %SF_check = EM2SE(SF_s_EM(:,ii),times_s_EM{ii}(end),mu_EM,mu_SE,phi_EM_0,L_EM,L_SE,T_EM,T_SE);
    % if (abs(SF_check(1) - mu1_SE) < 10^-4) && (SF_check(2) > 3e-3)
        for jj = 1:size(states_s_EM{ii},2)% size(states_s_EM{ii},2) returns the number of states per cell
             states_s_EM_SE{ii}(:,jj) = EM2SE(states_s_EM{ii}(:,jj),times_s_EM{ii}(jj),mu_EM,mu_SE,phi_EM_0,L_EM,L_SE,T_EM,T_SE);
             times_s_EM_SE{ii}(:,jj) = times_s_EM{ii}(jj)*T_EM/T_SE;
        end
        SF_s_EM_SE(:,ii) = states_s_EM_SE{ii}(:,end);
        %nn_EM_s_plus = nn_EM_s_plus+1;
       % end
        ii
end
    
% %Plots    
% figure    
% for ii = 1:size(states_s_EM,2)-1%nn_EM_s_plus-1
%         %figure
%         SCline = plot(states_s_EM_SE{ii}(1,:),states_s_EM_SE{ii}(2,:),'b');  SCline.Color(4) = 0.3;
%         hold on
% end
%     
% states_Moon_0 = EM2SE([params_EM.mu1;0;0;0],0, params_EM.mu_EM, params_EM.mu_SE, params_EM.phi_EM_0, params_EM.L_EM, params_EM.L_SE, params_EM.T_EM, params_EM.T_SE);
% moon_orbit_SE = circle(mu1_SE,0,L_EM/L_SE);
%     
% %figure
% grid on
% axis equal
% plot(mu1_SE,0,'ko','markerfacecolor','k','markersize',6)%interp2(X,Y,U,-mu2,0),'ro')
% text(mu1_SE*1.0001,0,'$m_E$','interpreter','latex','fontsize',15,'fontweight','bold')
% plot(states_Moon_0(1),states_Moon_0(2),'ko','markerfacecolor','k','markersize',6)
% text(states_Moon_0(1)*0.9993,states_Moon_0(2),'$m_{M,0}$','interpreter','latex','fontsize',15,'fontweight','bold')
% title(strcat(['Earth-Moon stable manifold, $\varphi_0$ = ' num2str(phi_EM_0*180/pi) '$^\circ$']),'interpreter','latex','fontsize',15)
% plot(moon_orbit_SE(:,1),moon_orbit_SE(:,2),'k--','linewidth',1)
% xlabel('$\bar{x}$','interpreter','latex')
% ylabel('$\bar{y}$','interpreter','latex')
% ylim = get(gca,'ylim');
% plot([mu1_SE,mu1_SE],[0 ylim(2)],'k','Linewidth',2)
% text(mu1_SE*1.0001,ylim(2)*0.5,'$S_0$','interpreter','latex','fontsize',15,'fontweight','bold')
% %set(gca,'ticklabelinterpreter','latex','fontsize',15)
% % savefig(strcat(['EM_transfer_manifolds_phi0_' num2str(phi_EM_0*180/pi) '.fig']))
    
