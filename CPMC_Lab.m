function [E_ave,E_err,Sz_ave, Sz_err,savedFileName]=CPMC_Lab(Lx,Ly,Lz,N_up,N_dn,kx,ky,kz,U,U_imp,tx,ty,tz,tx2,ty2,tz2,deltau,N_wlk,N_blksteps,N_eqblk,N_blk,itv_modsvd,itv_pc,itv_Em,t_bp,t_pop,suffix)
% function [E_ave,E_err,savedFileName]=CPMC_Lab(Lx,Ly,Lz,N_up,N_dn,kx,ky,kz,U,tx,ty,tz,deltau,N_wlk,N_blksteps,N_eqblk,N_blk,itv_modsvd,itv_pc,itv_Em, suffix)
% Perform a constrained path Monte Carlo calculatiion. Main function in the CPMC-Lab package
% Input
%   Lx: The number of lattice sites in the x direction.
%   Ly: The number of lattice sites in the y direction.
%   Lz: The number of lattice sites in the z direction.
%   N_up: The number of spin-up electrons
%   N_dn: The number of spin-down electrons
%   kx: The x component of the twist angle in TABC (twist-averaging boundary condition)
%   ky: The y component of the twist angle in TABC
%   kz: The z component of the twist angle in TABC
%   U: The on-site repulsion strength in the Hubbard Hamiltonian
%   tx: The hopping amplitude between nearest-neighbor sites in the x direction
%   ty: The hopping amplitude between nearest neighbor sites in the y direction
%   tz: The hopping amplitude between nearest neighbor sites in the z direction
%   deltau: The imaginary time step
%   N_wlk: The number of random walkers
%   N_blksteps: The number of random walk steps in each block
%   N_eqblk: The number of blocks used to equilibrate the random walk before energy measurement takes place
%   N_blk: The number of blocks used in the measurement phase
%   itv_modsvd: The interval between two adjacent modified Gram-Schmidt re-orthonormalization of the random walkers.
%   itv_pc: The interval between two adjacent population controls
%   itv_Em: The interval between two adjacent energy measurements
%   suffix: an identifying string e.g. timestamp to be appended to the end of the saved *.mat file
% Output:
%   E_ave: the ground state energy
%   E_err: the standard error in the ground state energy
%   savedFileName: name of the saved data file
%
% Huy Nguyen, Hao Shi, Jie Xu and Shiwei Zhang
% ©2014 v1.0
% Package homepage: http://cpmc-lab.wm.edu
% Distributed under the <a href="matlab: web('http://cpc.cs.qub.ac.uk/licence/licence.html')">Computer Physics Communications Non-Profit Use License</a>
% Any publications resulting from either applying or building on the present package 
%   should cite the following journal article (in addition to the relevant literature on the method):
% "CPMC-Lab: A Matlab Package for Constrained Path Monte Carlo Calculations" Comput. Phys. Commun. (2014)

%% Initialization
tic; % start the  timer
initialization; % initialize internal constants, form the trial wave function and assemble the initial population of walkers
format long;

N=floor(t_bp/N_blksteps);
flag_mea=0; %determine when a measurement should take place
flag_begin_bp=0; %determine when a back propagated measurement take place
flag_bp=0; %determine when the back propagating time ends
flag_bbp=zeros(N+1,1); %determine if back propagation is being effectuated 
E=0;
W=0;
Sz=zeros(N_sites,1);
% Preallocate arrays:
E_blk=zeros(N_blk,1); % array to store the energy measured in every block
W_blk=zeros(N_blk,1); % array to store the total weight in every block
% Preallocate back propagation-related arrays
Phi_old=zeros(N_sites, N_par, N_wlk, N+1);
parents=zeros(N_wlk,N+1);
B_up=zeros(t_bp+1, N_wlk, N_sites, N+1);
B_dn=zeros(t_bp+1, N_wlk, N_sites, N+1);
% Estimate how many back propagations are possible for the total proyection time
S_z_blk=zeros(N_blk-N-1,N_sites);
W_bp=zeros(N_blk-N-1,1);
% Measurement related parameters
ind_measure=0;
ind_parents=0;
ind_end=0;
list=zeros(N+1,1);
%% Equilibration phase
for i_blk=1:N_eqblk
    for j_step=1:N_blksteps
        [Phi, Phi_old, B_up, B_dn, parents, w, O, Sz, E, W, W_bp] = stepwlk(Phi, Phi_old, B_up, B_dn, parents, N_wlk, N_sites, N, w, O, Sz, E, W, W_bp, H_k, Proj_k_half_up, Proj_k_half_dn, flag_mea, flag_begin_bp, flag_bp, flag_bbp, Phi_T, N_up, N_par, t_bp, t_pop, U, U_imp, Lx, Ly, fac_norm, aux_fld, ind_parents, ind_end, list);
        if mod(j_step,itv_modsvd)==0
            [Phi, O] = stblz(Phi, N_wlk, O, N_up, N_par); % re-orthonormalize the walkers
        end
        if mod(j_step,itv_pc)==0
            [Phi, w, O, B_up, B_dn, parents]=pop_cntrl(Phi, w, O, B_up, B_dn, parents, N, N_wlk, N_sites, N_par, t_bp, flag_bbp); % population control
        end
    end
end

%% Measurement phase 
for i_blk=1:N_blk
    for j_step=1:N_blksteps
        t=(i_blk-1)*N_blksteps+j_step;
        if mod(j_step,itv_Em)==0
            flag_mea=1;
            if i_blk<N_blk-N
                flag_begin_bp=1;
                ind_parents=mod(ind_parents+1,N+1);
                flag_bbp(ind_parents+1)=1;
            end
        else
            flag_mea=0;
            flag_begin_bp=0;
        end
        if t>(N+1)*N_blksteps
            if j_step==t_bp-N*N_blksteps
                flag_bp=1;
                ind_measure=ind_measure+1;
                ind_end=mod(ind_end+1,N+1);
            else
                flag_bp=0;
            end
        end     
        for ii=1:N+1
            if flag_bbp(ii)==1
                list(ii)=list(ii)+1;
            else
                list(ii)=0;
            end
        end
        % propagate the walker
        if ind_measure>0
            [Phi, Phi_old, B_up, B_dn, parents, w, O, S_z_blk(ind_measure,:), E_blk(i_blk), W_blk(i_blk), W_bp(ind_measure)] = stepwlk(Phi, Phi_old, B_up, B_dn, parents, N_wlk, N_sites, N, w, O, S_z_blk(ind_measure,:), E_blk(i_blk), W_blk(i_blk), W_bp(ind_measure), H_k, Proj_k_half_up, Proj_k_half_dn, flag_mea, flag_begin_bp, flag_bp, flag_bbp, Phi_T, N_up, N_par, t_bp, t_pop, U, U_imp, Lx, Ly, fac_norm, aux_fld, ind_parents+1, ind_end+1, list);
        else
            [Phi, Phi_old, B_up, B_dn, parents, w, O, Sz, E_blk(i_blk), W_blk(i_blk), W] = stepwlk(Phi, Phi_old, B_up, B_dn, parents, N_wlk, N_sites, N, w, O, Sz, E_blk(i_blk), W_blk(i_blk), W, H_k, Proj_k_half_up, Proj_k_half_dn, flag_mea, flag_begin_bp, flag_bp, flag_bbp, Phi_T, N_up, N_par, t_bp, t_pop, U, U_imp, Lx, Ly, fac_norm, aux_fld, ind_parents+1, ind_end+1, list);
        end
        if mod(j_step,itv_modsvd)==0
            [Phi, O] = stblz(Phi, N_wlk, O, N_up, N_par); % re-orthonormalize the walkers
        end
        if mod(j_step,itv_pc)==0
            [Phi, w, O, B_up, B_dn, parents]=pop_cntrl(Phi, w, O, B_up, B_dn, parents, N, N_wlk, N_sites, N_par, t_bp, flag_bbp); % population control
        end
        if mod(j_step, itv_Em)==0
            % update the exponent of the pre-factor exp(-deltau*(H-E_T))
            fac_norm=(real(E_blk(i_blk)/W_blk(i_blk))-0.5*U*N_par)*deltau;
        end
        if flag_bp==1
            flag_bbp(ind_end+1)=0;    
        end
    end
    E_blk(i_blk)=E_blk(i_blk)/W_blk(i_blk);
    if ind_measure>0        
       S_z_blk(ind_measure,:)=S_z_blk(ind_measure,:)/W_bp(ind_measure);
    end
    display(strcat('E(',int2str(i_blk),')=',num2str(real(E_blk(i_blk)))))
end

%% Results
E=real(E_blk);
E_ave=mean(E)
E_err=std(E)/sqrt(N_blk)
Sz=real(S_z_blk);
Sz_ave=mean(Sz)
Sz_err=std(Sz)/sqrt(N_blk)
% The total computational time:
time=toc() % stops the timer
%% Save data to a *.mat file
save (savedFileName, 'E', 'E_ave', 'E_err', 'time');
save (savedFileName, '-append', 'Lx', 'Ly','Lz', 'N_up', 'N_dn', 'kx', 'ky','kz', 'U', 'tx', 'ty','tz');
save (savedFileName, '-append', 'deltau', 'N_wlk', 'N_blksteps', 'N_eqblk', 'N_blk', 'itv_pc','itv_modsvd','itv_Em');
save (savedFileName, '-append', 'H_k', 'Phi_T');

%% Explanation of saved quantities:
% E: the array of energy of each block
% time: The total computational time
% E_nonint_v: the non-interacting energy levels of the system
% Phi_T: the trial wave function
% For other saved quantities, type "help CPMC_Lab"
end