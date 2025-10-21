%% Main Code for Continuous-Time Version of Payment-Chain Model
% Uses Exact Solution
% (c) Bigio et al
%
% October 2025
% Main block happens here
% written by S. Bigio
clear; close all;

%% Running Preferences
plotiter=0;

%% Parameters (create subroutine)
% Model Parameters - Preferences and Technology
CT_parameters;

% Threshold values for q that lead to convergence
N_q=500;
q_vec=linspace(1,2,N_q);
b_h_mat=ones(N_q,1);

%% Analytic Solution
% solve consumption at deleveraging phase:
CT_solution;

% CT analytic Plots
CT_analyticplots;

%% Iterate over values of q such that no convergence takes place
for qq=1:N_q
    q=q_vec(qq);
    CT_solution;
    b_h_mat(qq)=b_h;
end

figure('Name','Hysteresis Threshold (Partial Equilibrium)')
plot(q_vec, b_h_mat,'LineWidth',2);
grid on; axis tight;
