%% Main Code for Continuous-Time Version of Payment-Chain Model
% Uses Exact Solution
% (c) Bigio et al
%
% October 2025
%%%%%%%%%%%%%%%%%%%%%%%
% computes continuous-time solution and threshold value that attracts
% convergence
%%%%%%%%%%%%%%%%%%%%%%%
% written by S. Bigio
clear; close all;

%% Running Preferences
plotiter=0;

%% Parameters (create subroutine)
% Grid values values for delta and q
N_delta=100;
N_q=500; 
q_max=12;
q_vec=linspace(1,q_max,N_q);
delta_vec=linspace(0.16,0.96,N_delta);

% Model Parameters - Preferences and Technology
CT_parameters;

% Initiating matrix
b_h_mat=zeros(N_delta,N_q);
b_h_delta=zeros(N_delta,1);

%% Main loop
for dd=1:N_delta
    delta=delta_vec(dd);
    for qq=1:N_q
        q=q_vec(qq);
        CT_solution;
        b_h_mat(dd,qq)=b_h;
    end
    mu_vec=1+rho*b_h_mat(dd,:);
    q_out=(CalA(mu_vec,delta)).^(-1);
    [~,index_delta]=min(abs(q_vec-q_out));
    b_h_delta(dd)=b_h_mat(dd,index_delta);
end
figure('Name','Hysteresis Threshold (Partial Equilibrium)')
plot(delta_vec, -b_h_delta,'LineWidth',2);
grid on; axis tight;

% Save thresholds from CT solution
delta_vec_ct=delta_vec; b_h_delta_ct=b_h_delta; 
save bh_cutoffs_ct.mat delta_vec_ct b_h_delta_ct;

function A_out=CalA(mu,delta) 
    if mu==1
        A_out=10e-16;
    elseif mu==0
        A_out=delta;
    else
        A_out=delta/(1-delta)*(1-mu)./mu.*log((1-mu.*delta)./(1-mu));
    end
end