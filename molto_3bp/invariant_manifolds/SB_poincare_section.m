      
nn_SE_u_plus = 1;
exit_phiB = 0;
    for ii = 1:size(states_u_SE,2)
        if states_u_SE{ii}(1,end) < pos_SE(2,1)% Stable manifold EM system (yplus part)
            for jj = 1:size(states_u_SE{ii},2)
               phi_B_t = atan2(states_u_SE{ii}(2,jj),(states_u_SE{ii}(1,jj)-params_SE.mu1));
                
                if (phi_B_t < phi_B) || states_u_SE{ii}(2,jj) < 0
                    states_u_SE_phiB{nn_SE_u_plus}(:,jj) = states_u_SE{ii}(:,jj);
                    times_u_SE_phiB{nn_SE_u_plus}(jj) = times_u_SE{ii}(jj);
                    exit_phiB = 1;
                end
            end
            SF_u_SE_phiB(:,nn_SE_u_plus) = states_u_SE_phiB{nn_SE_u_plus}(:,end);
            nn_SE_u_plus = nn_SE_u_plus+1;
%         else
%         states_u_SE_inf{nn_SE_u_inf} = states_u_SE{ii};
%         nn_SE_u_inf = nn_SE_u_inf+1; 
        end
        ii
    end

%Plots
% Plot1:
figure
    for ii = 1:nn_SE_u_plus-1
        %figure(3)
        SAline = plot(states_u_SE_phiB{ii}(1,:),states_u_SE_phiB{ii}(2,:),'r'); SAline.Color(4) = 0.3;
        hold on
    end
    
    %figure(3)
    plot(states_po_SE(1,:),states_po_SE(2,:),'k','Linewidth',2)
    grid on
    plot(mu1_SE,0,'ko','markerfacecolor','k','markersize',4)%interp2(X,Y,U,-mu2,0),'ro')
    text(mu1_SE*1.0001,0,'$m_E$','interpreter','latex','fontsize',15,'fontweight','bold')
    plot(pos_SE(2,1),pos_SE(2,2),'ko','markerfacecolor','k','markersize',4)
    text(pos_SE(2,1)*1.0001,pos_SE(2,2),'$L_2$','interpreter','latex','fontsize',15,'fontweight','bold')
    title(strcat(['Sun-Earth unstable manifold, $\varphi_B$ = ' num2str(phi_B*180/pi) '$^\circ$']),'interpreter','latex','fontsize',15)
    xlabel('$\bar{x}$','interpreter','latex')
    xlabel('$\bar{y}$','interpreter','latex')
    ang([mu1_SE,0],2e-3,[0 phi_B])
    
    if phi_B*180/pi>=90
        plot([mu1_SE,min(SF_u_SE_phiB(1,:))],[0, max(SF_u_SE_phiB(2,:))],'k','linewidth',1.5)
    else
        plot([mu1_SE,max(SF_u_SE_phiB(1,:))],[0, min(SF_u_SE_phiB(2,:))],'k','linewidth',1.5)
    end
    
    text(mu1_SE+2e-3,1e-3,'$\varphi_B$','interpreter','latex','fontsize',15,'fontweight','bold')
    %set(gca,'ticklabelinterpreter','latex','fontsize',15)
    %savefig(strcat(['SB_manifolds' num2str(phi_B*180/pi) '.fig']))
    

%%%%Para qué es estoo??? y dónde está Moon capture leg?? 
%%Plotea phase space de la pagina 39?

% nn_SE_u_plus = 1;
% nn_SE_u_min = 1;
% for ii = 1:size(states_s_SE,2)
%     if SF_u_SE(2,ii) > 0
%         SF_u_SE_yplus(:,nn_SE_u_plus) = SF_u_SE(:,ii);
%         nn_SE_u_plus = nn_SE_u_plus+1;
%     else
%         SF_u_SE_ymin(:,nn_SE_u_min) = SF_u_SE(:,ii);
%         nn_SE_u_min = nn_SE_u_min+1;
%     end
% end
% 
% figure
% %     y_s_EM_yplus_SE = SF_s_EM_yplus_SE(2,:)';
% %     ydot_s_EM_yplus_SE = SF_s_EM_yplus_SE(4,:)';
% %     k = boundary(y_s_EM_yplus_SE,ydot_s_EM_yplus_SE,0.1);
% %     plot(y_s_EM_yplus_SE,ydot_s_EM_yplus_SE,'o','color','b','markersize',2,'markerfacecolor','b')
% %     hold on
% %     plot(y_s_EM_yplus_SE(k),ydot_s_EM_yplus_SE(k),'color','b')
% y_u_SE_yplus = SF_u_SE_yplus(2,:)';
% ydot_u_SE_yplus = SF_u_SE_yplus(4,:)';
% k = boundary(y_u_SE_yplus,ydot_u_SE_yplus,0.85);
% plot(y_u_SE_yplus,ydot_u_SE_yplus,'o','color','r','markersize',2,'markerfacecolor','b')
% hold on
% plot(y_u_SE_yplus(k),ydot_u_SE_yplus(k),'color','r')
% hold on
% y_s_SE_ymin = SF_s_SE_ymin(2,:)';
% ydot_s_SE_ymin = SF_s_SE_ymin(4,:)';
% k = boundary(y_s_SE_ymin,ydot_s_SE_ymin,0.85);
% plot(y_s_SE_ymin,ydot_s_SE_ymin,'o','color','b','markersize',2,'markerfacecolor','b')
% hold on
% plot(y_s_SE_ymin(k),ydot_s_SE_ymin(k),'color','b')
% hold on
% 
% xlabel('$\bar{y}$','interpreter','latex')
% ylabel('$\dot{y}$','interpreter','latex')
% title('Phase space in Poincare section','interpreter','latex')
% grid on
% %savefig('SE_phspc_patch.fig')
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% nn_SE_s_plus = 1;
% nn_SE_s_min = 1;
% 
% for ii = 1:size(states_s_SE,2)
%     if SF_s_SE(2,ii) > 0
%         % Stable manifold EM system (yplus part)
%         SF_s_SE_yplus(:,nn_SE_s_plus) = SF_s_SE(:,ii);
%         states_s_EM{ii}(1,:)
%         
%         nn_SE_s_plus = nn_SE_s_plus+1;
%     else
%         % Stable manifold EM system (yminus part)
%         SF_s_SE_ymin(:,nn_SE_s_min) = SF_s_SE(:,ii);
%         nn_SE_s_min = nn_SE_s_min+1;
%     end
% end
