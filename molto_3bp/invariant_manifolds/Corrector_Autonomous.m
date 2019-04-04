function [X0,T0,Error,Floquet]=Corrector_Autonomous(Fun_DF,X0,T0,Itmax,Tol,TolRel,TolAbs,dh,Ind_Fix)

% Periodic Orbit Corrector for autonomous system
% Fun_DF                    - > function F = Fun(t,x) that gives the right-hand side 
% X0                        - > Initial Condition guess
% T0                        - > Period guess
% Itmax, Tol, TolRel,TolAbs - > Numerical Parameters
% dh                        - > Jacobian: Finite differences step 
% Ind_Fix                   - > (Integer) Index of the state vector variable that is fixed

%%%%%%%%%%%%%%%
%% Example %%%%
%%%%%%%%%%%%%%%
% Let the system 
% x' = x*y
% y' = 2x-y^2
% If you make the call
% [X0,T0,Error,Floquet]=Corrector_Autonomous(Fun_DF,X0,T0,Itmax,Tol,TolRel,TolAbs, 2)
% The program will try to find a periodic orbit that at t=0 satisfies y(t=0) = X0(2)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N        = length(X0);
Iden     = eye(N);
options  = odeset('RelTol',TolRel,'AbsTol',ones(1,N*(1+N))*TolAbs);
options2 = odeset('RelTol',TolRel,'AbsTol',ones(1,N)*TolAbs);

XFix     = X0(Ind_Fix); 

for i=1:1:Itmax
  
    % Compute the Initial Error
    X0(Ind_Fix) = XFix; 
    [T X]       = ode45(Fun_DF,[0 T0],X0,options2);
    Error       = max(abs(X(1,:)-X(end,:)));
    display(['Corrector: iteracion = ' num2str(i) '  Error = ' num2str(Error)])
         
    
    Xvar = [];
    Xvar(1,1:N) = X0';
    for j=1:1:N
        Xvar(1,N+(j-1)*N+1:N+j*N) = Iden(j,:);
    end
    Xvar = Xvar';
        
    if Error<Tol      
        
        if i==1       
           [T X]    = ode45(@FUNvariationalDF,[0 T0],Xvar,options);
            for j=1:1:N
               for k=1:1:N
                M(j,k)   = X(end,N+(j-1)*N+k);
               end
            end
            M=M';
            [Vec Val] = eig(M);
            for j=1:1:N
                Floquet(1,j)  = Val(j,j);
            end  
        end
        break
    end
   
    % Compute Monodromy Matrix
    [T X]     = ode45(@FUNvariationalDF,[0 T0],Xvar,options);
    XF        = X(end,1:N*(N+1))';
    for j=1:1:N
        for k=1:1:N
            M(j,k) = X(end,N+(j-1)*N+k);
        end
    end
    M=M';
    % Compute the derivative 
    DF = feval(Fun_DF,T0,X0);
    
    % Prepare equation A * [DeltT X1 X2... XN ] = b (without X(Index))
    A            = M-eye(N);
    A(:,Ind_Fix) = DF;
    B            = -(XF(1:N)-X0);
    Correc       = A\B;
    for j=1:1:N
        if j==Ind_Fix
            T0 = T0+Correc(j);
        else
            X0(j,1)     = X0(j,1)+Correc(j);
        end
    end
    [Vec Val] = eig(M);
    for j=1:1:N
        Floquet(1,j)   = Val(j,j);
    end
    if i==Itmax
         [T X]     = ode45(Fun_DF,[0 T0],X0,options2);
         Error = max(abs(X(1,:)-X(end,:)));
         display(['Corrector: iteracion = ' num2str(i) '  Error = ' num2str(Error)])
    end
end



  function  DF = FUNvariationalDF(t,YV)
        
        Jac       = Fun_Jac_Num(Fun_DF,t,YV(1:N,1),dh);
        DF(1:N,1) = feval(Fun_DF,t,YV(1:N,1));   

          
        for cont1 = 1:1:N                                                            % Loop for Xi variational
            for cont2 = 1:1:N                                                        % Index of the Xi
               DF(N + (cont1-1)*N + cont2) = 0;
               for cont3 = 1:1:N                                                    % Loop for the sum
                 DF(N + (cont1-1)*N + cont2) = DF(N + (cont1-1)*N + cont2) + Jac(cont2,cont3)*YV(N + N*(cont1-1) + cont3);
               end
            end
        end
    end

end