% Get Poincaré section on S_A with varying phi_A (SE stable manifold)
             
nn_SE_s_min = 1;
exit_phiA = 0;
for ii = 1:size(states_s_SE,2)
  if states_s_SE{ii}(1,end) < pos_SE(2,1)% Stable manifold EM system (yplus part)
     for jj = 1:size(states_s_SE{ii},2)
      phi_A_t = atan2(-states_s_SE{ii}(2,jj),(states_s_SE{ii}(1,jj)-params_SE.mu1));
                    
        if (phi_A_t < phi_A) || states_s_SE{ii}(2,jj) > 0
          states_s_SE_phiA{nn_SE_s_min}(:,jj) = states_s_SE{ii}(:,jj);
          times_s_SE_phiA{nn_SE_s_min}(jj) = times_s_SE{ii}(jj);
          exit_phiA = 1;
          % elseif phi_A_t > phi_A && exit_phiA == 1
          % SF_s_SE_phiA(:,nn_SE_s_min) = states_s_SE_phiA{nn_SE_s_min}(:,end);
        end
      end
          SF_s_SE_phiA(:,nn_SE_s_min) = states_s_SE_phiA{nn_SE_s_min}(:,end);
          nn_SE_s_min = nn_SE_s_min+1;
%     elseif states_s_SE{ii}(1,end) > pos_SE(2,1)
%         states_s_SE_inf{nn_SE_s_inf} = states_s_SE{ii};
%         nn_SE_s_inf = nn_SE_s_inf+1;            
  end
   ii
end

% save(strcat(['states_s_SE_phiA_' num2str(phi_A*180/pi)]),'states_s_SE_phiA')
% save(strcat(['times_s_SE_phiA_' num2str(phi_A*180/pi)]),'times_s_SE_phiA')
% save(strcat(['SF_s_SE_phiA_' num2str(phi_A*180/pi)]),'SF_s_SE_phiA')
   
%Plots
  figure      
 for ii = 1:nn_SE_s_min-1
    %figure(3)
    SAline = plot(states_s_SE_phiA{ii}(1,:),states_s_SE_phiA{ii}(2,:),'b'); SAline.Color(4) = 0.3;
    hold on
 end
        
 %figure(3)
 plot(states_po_SE(1,:),states_po_SE(2,:),'k','Linewidth',2)
 grid on
 plot(mu1_SE,0,'ko','markerfacecolor','k','markersize',4)%interp2(X,Y,U,-mu2,0),'ro')
 text(mu1_SE*1.0001,0,'$m_E$','interpreter','latex','fontsize',15,'fontweight','bold')
 plot(pos_SE(2,1),pos_SE(2,2),'ko','markerfacecolor','k','markersize',4)
 text(pos_SE(2,1)*1.0001,pos_SE(2,2),'$L_2$','interpreter','latex','fontsize',15,'fontweight','bold')
 title(strcat(['Sun-Earth stable manifold, $\varphi_A$ = ' num2str(phi_A*180/pi) '$^\circ$']),'interpreter','latex','fontsize',15)
 xlabel('$\bar{x}$','interpreter','latex')
 xlabel('$\bar{y}$','interpreter','latex')
 ang([mu1_SE,0],2e-3,[2*pi-phi_A 2*pi])
 
 if phi_A*180/pi>90
     plot([mu1_SE,min(SF_s_SE_phiA(1,:))],[0, min(SF_s_SE_phiA(2,:))],'k','linewidth',1.5)
 else
     plot([mu1_SE,max(SF_s_SE_phiA(1,:))],[0, min(SF_s_SE_phiA(2,:))],'k','linewidth',1.5)
 end
 text(mu1_SE+2e-3,-1e-3,'$\varphi_A$','interpreter','latex','fontsize',15,'fontweight','bold')
 %set(gca,'ticklabelinterpreter','latex','fontsize',15)
 %savefig(strcat(['SA_manifolds' num2str(phi_A*180/pi) '.fig']))
        
 

    
