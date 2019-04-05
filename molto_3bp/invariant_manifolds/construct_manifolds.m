function [states_s, times_s, SF_s, states_u, times_u, SF_u] = construct_manifolds(params, T0, states_po, times_po, eigvec, eigval, inv_phi_0, prnt_out_dt, npoints, stop_fun)


%% Step 4: Computation of invariant manifolds
T_po = T0;
idx = round(linspace(1,numel(times_po),npoints));
et0 = 0;
eps = 10^-4;
deltat = 5*T_po;%why 5?
sign = [-1 1];

for ii = 1:npoints
    for jj =1:numel(sign)
    X = states_po(:,idx(ii));
    t = times_po(idx(ii));
    phi = [eigvec(:,1)*exp(eigval(1)*t) eigvec(:,2)*exp(eigval(2)*t) eigvec(:,3)*exp(eigval(3)*t) eigvec(:,4)*exp(eigval(4)*t)];
    Phi = phi*inv_phi_0;
    [mon_eigv, mon_eigs] = eig(Phi);

    Yu = mon_eigv(:,1); % unstable eigenvector
    Ys = mon_eigv(:,2); % stable eigenvector

    Xu = real(states_po(:,idx(ii)) + eps*sign(jj)*Yu);
    Xs = real(states_po(:,idx(ii)) + eps*sign(jj)*Ys);
    
    %Integrate unstable manifold forwards in time
    [SF, etf_u, states, times] = PCR3BP_propagator (Xu, params,...
        et0, deltat, prnt_out_dt, stop_fun);
    states_u{2*(ii-1)+jj} = states;
    times_u{2*(ii-1)+jj} = times;
    SF_u(:,2*(ii-1)+jj) = SF;
    
    
    %Integrate stable manifold backwards in time    
   [SF, etf, states, times] = PCR3BP_propagator (Xs, params,...
       et0, -deltat, prnt_out_dt, stop_fun);
   states_s{2*(ii-1)+jj} = states;
   times_s{2*(ii-1)+jj} = times;
   SF_s(:,2*(ii-1)+jj) = SF;

    2*(ii-1)+jj%display iteration
    end
end
end