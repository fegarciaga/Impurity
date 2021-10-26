function [spin_z]=measure_bp(Phi_T, phi_old, B_up, B_dn, Proj_k_half_up, Proj_k_half_dn, t_bp, t_pop, N_up, N_par, N_sites)
  %% Implement backwards propagation with the stored values
  % calculate backpropagated bra
  Phi_left=Phi_T;
    for ii=1:t_bp+1
        Phi_left(:,1:N_up)=Proj_k_half_up*Phi_left(:,1:N_up);
        Phi_left(:,1+N_up:N_par)=Proj_k_half_dn*Phi_left(:,1+N_up:N_par);
        for jj=1:N_sites
            Phi_left(jj,1:N_up)=Phi_left(jj,1:N_up)*B_up(t_bp-ii+2,jj);
            Phi_left(jj,N_up+1:N_par)=Phi_left(jj,1+N_up:N_par)*B_dn(t_bp-ii+2,jj);
        end
        Phi_left(:,1:N_up)=Proj_k_half_up*Phi_left(:,1:N_up);
        Phi_left(:,1+N_up:N_par)=Proj_k_half_dn*Phi_left(:,1+N_up:N_par);
        if mod(ii,t_pop)==0
            [Phi_left(:,1:N_up), Trsh]=qr(Phi_left(:,1:N_up),0);
            [Phi_left(:,1+N_up:N_par), Trsh]=qr(Phi_left(:,1+N_up:N_par),0);
        end
    end
    %% Once the bra is back propagated, the observable is calculated the usual way
    % calculate inverse matrices in order to calculate Green's function
    %% Check when it doesn't work 
    inv_O_up=inv(Phi_left(:,1:N_up)'*phi_old(:,1:N_up));
    inv_O_dn=inv(Phi_left(:,1+N_up:N_par)'*phi_old(:,1+N_up:N_par));
    temp_up=phi_old(:,1:N_up)*inv_O_up;
    temp_dn=phi_old(:,N_up+1:N_par)*inv_O_dn;
    G_up=temp_up*Phi_left(:,1:N_up)';
    G_dn=temp_dn*Phi_left(:,N_up+1:N_par)';
    Id=eye(N_sites);
    G_up_tilde=Id-G_up;
    G_dn_tilde=Id-G_dn;
    % As second try I will measure the spin correlation
    aux1=0;
    aux2=0;
    for i=1:N_sites
        % Measures the spin for each site
        spin_z(i)=0.5*(G_up(i,i)-G_dn(i,i));
    end
end