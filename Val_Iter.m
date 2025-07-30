%% Payment Chains
% (c) 
% version with only consumption, no production. 
clear; close all;
params.tol=10e-6;
params.printinfo=0;
params.printit=0;

%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%
% Value functions - Partial Equilibrium
epsilon=10e-10;
params.N_b=1000;
Val_Iter_params;

% Support of values of debt
range=[-1 params.b_bar-epsilon]    ;
B_vec=linspace(range(1),range(end),params.N_b); % grid space
params.B_vec=B_vec;
params.range=range;

%% Load Functions
Val_Iter_functions;

%% Main Iteration
[V0,policies,B_exit_num]=solveValueFunction(B_vec,funcs,params);

%% Computing Conjectured Value after which you do not exit
% Tv_iter=Tv(B_vec,Bast);
% [V_Bast,I]=max(Tv_iter(Tv_iter==real(Tv_iter)));

%% Main Plots
Val_Iter_plots(params, funcs, V0, policies, B_exit_num);

%% Simulated Paths
T=25; % periods
init_bs = (350:100:850); % initial conditions
% init_bs = [10 N_b-10 N_b/2]; % initial conditions

% compute histories
[hist]=Val_iter_sim(T,init_bs,B_vec,policies);

% plot histories
Val_Iter_simplots(hist);

%% Val Iter Tests
Val_Iter_tests;

%% Functions
function [V_out,policies,B_exit_num]=solveValueFunction(B_vec,funcs,params)
    % unpack parameters
    printinfo=params.printinfo;
    B_tilde=params.B_tilde;
    tol=params.tol;
    E=funcs.E;
    C_w=funcs.C_w;
    S_w=funcs.S_w;
    X_w=funcs.X_w;
    val_low=funcs.val_low;
    
    % initialize values
    V0=val_low(B_vec); V_out=V0; C_out=0*V_out; S_out=0*V_out; X_out=0*V_out; Bp_out=0*V_out;
    
    % reset tolerance to enter loop
    val_gap=2*tol;

    % while loop iteration
    while val_gap>tol
        % print information
        if printinfo
            plot(B_vec,V0,'LineWidth',1,'LineStyle',':'); drawnow;
            disp(['Val at iter:' num2str(val_gap)]);
        end

        % print information
        bb=0;
        for B=B_vec
            bb=bb+1;
            % maximize value function based on cubic spline
            Tv=@(B_p,B) log(C_w(E(B_p,B),B,B_tilde))+params.beta*interp1(B_vec,V0,B_p,'cubic');
            % compute all possible values
            Tv_iter=Tv(B_vec,B);
            index=(imag(Tv_iter)==0);
            % evaluate at optimum to obtain solution
            [V_out(bb),I]=max(Tv_iter(index));
            Bp_aux=B_vec(index);
            % update solutions
            C_out(bb)=C_w(E(Bp_aux(I),B),B,B_tilde);
            S_out(bb)=S_w(E(Bp_aux(I),B),B,B_tilde);
            X_out(bb)=X_w(E(Bp_aux(I),B),B,B_tilde);
            Bp_out(bb)=Bp_aux(I);
        end
        val_gap=max(abs(V_out-V0));
        V0=V_out;
    end
    % Numerical Exit point
    B_exit_num=min(B_vec(abs(V0-val_low(B_vec))<10e-10));
    
    % recording policies
    policies.C_out=C_out;
    policies.S_out=S_out;
    policies.X_out=X_out;
    policies.Bp_out=Bp_out;
end

