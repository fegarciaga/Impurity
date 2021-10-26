function [Phi, w, O, B_up, B_dn, parents]=pop_cntrl(Phi, w, O, B_up, B_dn, parents, N, N_wlk, N_sites, N_par, t_bp, flag_bbp)
% function [Phi, w, O]=pop_cntrl(Phi, w, O, N_wlk, N_sites, N_par)
% Perform population control with a simple "combing" method
% Inputs:
%   Phi: the whole ensemble of walkers
%   w: array containing the weights of all the walkers
%   O: array containing the overlaps of all the walkers
%   N_wlk: the number of walkers
%   N_sites: the total number of lattice sites
%   N_par: the total number of electrons
% Outputs:
%   Phi: the new ensemble of walkers after population control
%   w: the new array of weights
%   O: the new array of overlaps
%
% Huy Nguyen, Hao Shi, Jie Xu and Shiwei Zhang
% ©2014 v1.0
% Package homepage: http://cpmc-lab.wm.edu
% Distributed under the <a href="matlab: web('http://cpc.cs.qub.ac.uk/licence/licence.html')">Computer Physics Communications Non-Profit Use License</a>
% Any publications resulting from either applying or building on the present package 
%   should cite the following journal article (in addition to the relevant literature on the method):
% "CPMC-Lab: A Matlab Package for Constrained Path Monte Carlo Calculations" Comput. Phys. Commun. (2014)

    %% Preparation
    % Create empty matrices that will hold the outputs
    new_Phi=zeros(N_sites, N_par, N_wlk); %in the end the number of walkers will still be N_wlk
    new_O=zeros(N_wlk,1);
    new_B_up=zeros(t_bp+1, N_wlk, N_sites,N+1);
    new_B_dn=zeros(t_bp+1, N_wlk, N_sites,N+1);
    new_parents=zeros(N_wlk,N+1);
    % scaling factor to bring the current total weight back to the original
    % level (=N_wlk)
    d=N_wlk/sum(w);
    % start the "comb" at a random position to avoid bias against the first walker
    sum_w=-rand;
    n_wlk=0;

    %% Apply the comb
    for i_wlk=1:N_wlk
        sum_w=sum_w+w(i_wlk)*d;
        n=ceil(sum_w);
        for j=(n_wlk+1):n
            new_Phi(:,:,j)=Phi(:,:,i_wlk);
            new_O(j)=O(i_wlk);
            % for copied determinants, the ancestral path is actualized aswell as the corresponding father
            for i=1:N+1
                if flag_bbp(i)==1
                    new_B_up(:,j,:,i)=B_up(:,i_wlk,:,i);
                    new_B_dn(:,j,:,i)=B_dn(:,i_wlk,:,i);
                    new_parents(j,i)=parents(i_wlk,i);
                end
            end
        end
        n_wlk=n;
    end

    %% Return the new population, weights and overlaps:
    Phi=new_Phi;
    O=new_O;
    % All new walkers have weights to 1 and the total weight = N_wlk
    w=ones(N_wlk,1);
    B_up=new_B_up;
    B_dn=new_B_dn;
    parents=new_parents;
end