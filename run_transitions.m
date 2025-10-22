%% Payment Crises: Code to Compute Equilibrium Solution 
% (c) Saki Bigio
clear; close all;
clear; close all; 
addpath('functions', 'plotting', 'tests');
params.folder='output/';

% Code Specs
printit=1; % Save plots into PDF
plotit =1; % 1 if plot needed

% Plot Specs
color_t1=[0.7 0.7 0.7];
color_t0=[0.0 0.0 0.6];
linethick=3;
color_arrow=[0 0 0];
fontsize=16;

%% Parameter Set
Trans_params; 

% Values for computing thresholds
N_delta=100;
delta_vec=linspace(0.16,0.96,N_delta);

%% Equilibrium Functions
Trans_function;

%% Test Plots
% internal plots to test model
if plotit==2
    PlotTest_Equilibrium;
end

%% Transitional Dynamics
% Transitional Dynamics plots
B_tilde=2; 
B_tilde_p=B_tilde; % Initial Values for B_tilde
B_star_aux=Bstar(B_tilde);
B_ast_aux=Bast(B_tilde_p);
N_t=20;

% Baseline Version plots
if plotit==1
    Plots_AggregateEuler;
end

% Extreme version of plots
if plotit==0
    params.delta=0.5;
    Trans_function; % Important, update functions
    Plots_AggregateEuler;
    Trans_params;
    Trans_function;
end

%% Evaluate Different Thresholds - Delta
B_ast_vec=zeros(N_delta,1);
for dd=1:N_delta
    params.delta=delta_vec(dd);
    Trans_function; 
    B_ast_vec(dd)=Bast(B_tilde_p);
end
figure('Name','Hysteresis Region')
plot(delta_vec,B_ast_vec,'Color','k','LineStyle',':'); grid on; hold on;
fill([delta_vec delta_vec(1)],[B_ast_vec; B_ast_vec(end)],[0.8 0.8 0.99]); grid on;
xlabel('$\delta$','FontSize', 16,'Interpreter','latex');
ylabel('$B^{\ast}$','FontSize', 16,'Interpreter','latex');
text(delta_vec(1),B_ast_vec(end),'$\downarrow$ Efficient Region ', 'HorizontalAlignment','Left','VerticalAlignment','top','Rotation',90,'interpreter', 'latex','FontSize',12);
load("data/bh_cutoffs_ct.mat")
hold on;
plot(delta_vec_ct,-b_h_delta_ct,'Color','k','LineStyle',':'); grid on;
fill([delta_vec_ct delta_vec_ct(end)],[-b_h_delta_ct; -b_h_delta_ct(1)],[0.9 0.9 0.9]); grid on;
text(delta_vec(end)-0.02,B_ast_vec(1)*0.9,'$\uparrow$ Hysteresis Region ', 'HorizontalAlignment','Left','VerticalAlignment','top','Rotation',90,'interpreter', 'latex','FontSize',12);
text(delta_vec(end)-0.02,B_ast_vec(end)*0.99,'$\uparrow$ No Symmetric Eq.', 'HorizontalAlignment','Left','VerticalAlignment','top','Rotation',90,'interpreter', 'latex','FontSize',12);
axis tight;
if printit==1
    orient landscape;
    ax = gca;
    exportgraphics(ax,[params.folder 'F_delta_domains.pdf'],'BackgroundColor','none');
end

% Rest parameters
Trans_params; % reset parameters
Trans_function;

%% Steady-State in 3-d
% 3-d plots
if plotit==1
    Plots_SteadyState_3d;
end

%% Resilience Plot
% Variables
varlist={'B_tilde_t','B_t','R_t','C_s_t','S_w_t','mu_t','q_t','X_w_t','Y_t'};
vartaglist={'$\tilde{B}_t$','$B_t$','$R_t$','$C^s_t$','$S^w_t$','$\mu_t$','$q_t$','$X^w_t$','$Y_t$'};
printlist={'fig_tildeB','fig_B','fig_R','fig_C','fig_S','fig_mu','fig_q','fig_X','fig_Y'};

% Main Values for Reslience plot
B0_base     = 1.00              ; % initial value
B0          = B0_base           ; % Initial value for simuluations
B_tilde_rho = 0.6              ; % mean reversal rate after recovery
B_tilde_ss  = 0.1*params.B_bar ; % steady-state smooth borrowing limit
B_tilde_init= 0.1*B_tilde_ss   ; 
B_tilde_pre = B_tilde_ss        ;
B_star      = Bstar(B_tilde_ss) ;

% Construct Shock
T_pre      = 1  ; index_pre=1:T_pre;
T_crunch   = 30 ; index_crunch=T_pre+1:T_pre+T_crunch;
T_rec      = 80; index_rec=T_pre+T_crunch+1:T_pre+T_crunch+T_rec+1; 
T_post     = 10 ; index_post = T_pre + T_crunch + T_rec + 1 : T_pre + T_crunch + T_rec + T_post;

% Example #1: Smooth Transition  
B_tilde_s   = [ones(1,T_pre)*B_tilde_ss zeros(1,T_crunch)+B_tilde_pre B_tilde_ss-B_tilde_rho.^(0:T_rec)*(B_tilde_ss-B_tilde_init) ones(1,T_post)*B_tilde_ss]';
sim_smooth=sim_mod(B_tilde_s,params,funcs,B0);

% run same plot with lower delta
params.delta=0.85;
Trans_function; % Important, update functions
B0          = B0_base;
sim_deltachange=sim_mod(B_tilde_s,params,funcs,B0);
if plotit==0
    for vv=1:numel(varlist)
        var_baseline=sim_smooth.(varlist{vv});
        var_higherd =sim_deltachange.(varlist{vv});
        vartag=vartaglist{vv};
        figure('Name','Resilience Plot');
       %  title(vartag,'Interpreter','latex','FontSize',fontsize); 
        hold on; grid on; 
        plot_comp2(index_pre,index_crunch,index_rec,index_post,var_baseline,var_higherd);
        if vv==1
            legend('High Delta','Lower Delta','Interpreter','Latex','FontSize',fontsize,'Location','SouthEast','Box','off');
             if 1==1
                ylimits=ylim;
                text(index_pre(end)+1,ylimits(1)+0.05,'Deleveraging Phase', 'HorizontalAlignment','Left','VerticalAlignment','top','Rotation',90,'interpreter', 'latex','FontSize',12);
                text(index_crunch(end)+1,ylimits(1)+0.05,'Shock Phase', 'HorizontalAlignment','Left','VerticalAlignment','top','Rotation',90,'interpreter', 'latex','FontSize',12);
            end
        end
        if printit==1
            orient landscape;
            ax = gca;
            exportgraphics(ax,[params.folder printlist{vv} '_resilience.pdf'],'BackgroundColor','none');
        end
    end
end

% Resetting Parameters
Trans_params;
Trans_function;

%% Planner vs. Market Solution Comparison (smooth and violent transitions)
% Variables
varlist={'B_tilde_t','B_t','R_t','C_s_t','S_w_t','mu_t','q_t','X_w_t','Y_t'};
vartaglist={'$\tilde{B}_t$','$B_t$','$R_t$','$C^s_t$','$S^w_t$','$\mu_t$','$q_t$','$X^w_t$','$Y_t$'};
printlist={'fig_tildeB','fig_B','fig_R','fig_C','fig_S','fig_mu','fig_q','fig_X','fig_Y'};

% Main Values for Reslience plot
B_tilde_rho = 0.95             ; % mean reversal rate after recovery
B_tilde_ss  = 0.1*params.B_bar; % steady-state smooth borrowing limit
B_tilde_init= 0.0*B_tilde_ss   ; 
B_tilde_pre = B_tilde_init     ;
B0          = Bstar(B_tilde_ss); % Initial value for simuluations
B_star      = Bstar(B_tilde_ss);

% Construct Shock
T_pre      = 1  ; index_pre=1:T_pre;
T_crunch   = 20 ; index_crunch=T_pre+1:T_pre+T_crunch;
T_rec      = 80 ; index_rec=T_pre+T_crunch+1:T_pre+T_crunch+T_rec+1; 
T_post     = 1  ; index_post = T_pre + T_crunch + T_rec + 1 : T_pre + T_crunch + T_rec + T_post;

% Example #1: Violent Transitions 
% Market - Violent Transition
B_tilde_v   = [ones(1,T_pre)*B_tilde_ss zeros(1,T_crunch)+B_tilde_init ones(1,+T_rec+T_post)*B_tilde_ss]';
sim_violent=sim_mod(B_tilde_v,params,funcs,B0);

% Planner - Violent Transition
B_ss_violent=sim_violent.B_t(end);
theta_violent=1-(1-params.beta)*B_ss_violent;
sim_p_violent=sim_planner(B_tilde_v,theta_violent,params,funcs,B0);

% Example #2: Smooth Transition 
% Market - Smooth Transition
B_tilde_s   = [ones(1,T_pre)*B_tilde_ss zeros(1,T_crunch)+B_tilde_pre B_tilde_ss-B_tilde_rho.^(1:T_rec+1)*(B_tilde_ss-B_tilde_init) ones(1,T_post)*B_tilde_ss]';
sim_smooth=sim_mod(B_tilde_s,params,funcs,B0);
% Planner  - Smooth Transition
B_ss_smooth=sim_smooth.B_t(end);
theta_smooth=1-(1-params.beta)*B_ss_smooth;
sim_p_smooth=sim_planner(B_tilde_s,theta_smooth,params,funcs,B0);

% Compare Violent and Smooth Transitions
if plotit==1
    for vv=1:numel(varlist)
        var_baseline=sim_smooth.(varlist{vv});
        var_higherd =sim_violent.(varlist{vv});
        vartag=vartaglist{vv};
        figure('Name','Comparison Plot');
        title(vartag,'Interpreter','latex','FontSize',fontsize); 
        hold on; grid on; 
        plot_comp2(index_pre,index_crunch,index_rec,index_post,var_baseline,var_higherd);
        if vv==1
            legend('Smooth','Violent','Interpreter','Latex','FontSize',fontsize,'Location','SouthEast','Box','off');
             if 1==1
                ylimits=ylim;
                text(index_pre(end)+1,ylimits(1)+0.05,'Deleveraging Phase', 'HorizontalAlignment','Left','VerticalAlignment','top','Rotation',90,'interpreter', 'latex','FontSize',12);
                text(index_crunch(end)+1,ylimits(1)+0.05,'Shock Phase', 'HorizontalAlignment','Left','VerticalAlignment','top','Rotation',90,'interpreter', 'latex','FontSize',12);
            end
        end
        if printit==1
            orient landscape;
            ax = gca;
            exportgraphics(ax,[params.folder printlist{vv} '.pdf'],'BackgroundColor','none');
        end
    end
end

% Main Simulation
printlist={'fig2_tildeB','fig2_B','fig2_R','fig2_C','fig2_S','fig2_mu','fig2_q','fig2_X','fig2_Y'};
if plotit==1
    for vv=1:numel(varlist)
        var_v_ce=sim_violent.(varlist{vv});
        var_v_rp=sim_p_violent.(varlist{vv});
        vartag=vartaglist{vv};
        figure('Name','Planner vs. Market (Violent)');
        hold on; grid on; 
        plot_comp2(index_pre,index_crunch,index_rec,index_post,var_v_ce,var_v_rp);
        title(vartag,'Interpreter','latex','FontSize',fontsize); 
        xlabel('Time','Interpreter','latex','FontSize',fontsize); title(vartag,'Interpreter','latex','FontSize',fontsize); hold on; 
        if vv==1
            legend('Violent Transition (CE)','Violent Transition (PE)','Interpreter','Latex','FontSize',fontsize,'Location','NorthWest','Box','off');
        end
        if printit==1
            orient landscape;
            eval(['saveas(gcf,''' printlist{vv} ''',''pdf'');']);
        end
    end
end


printlist={'fig3_tildeB','fig3_B','fig3_R','fig3_C','fig3_S','fig3_mu','fig3_q','fig3_X','fig3_Y'};
for vv=1:numel(varlist)
    var_s_ce=sim_smooth.(varlist{vv});
    var_s_rp=sim_p_smooth.(varlist{vv});
    vartag=vartaglist{vv};
    figure('Name','Comp-Planner (smooth shock)');
    hold on; grid on; 
   %  title(vartag,'Interpreter','latex','FontSize',fontsize); 
    plot_comp2(index_pre,index_crunch,index_rec,index_post,var_s_ce,var_s_rp);
    if vv==1
        legend('Competitive Equilbrium','Planner Solution','Interpreter','Latex','FontSize',fontsize,'Location','SouthEast','Box','off');
    end
    if printit==1
        orient landscape;
        ax = gca;
        exportgraphics(ax,[params.folder printlist{vv} '.pdf'],'BackgroundColor','none');
    end
end

%% Equilibrium Rate Function
function [R_out,fval,flag]=R(B,B_tilde,B_tilde_p,params)
    beta=params.beta;
    e_cond=(B_tilde>B&&B_tilde<1+beta*B); % expenditure condition
    s_cond=(B_tilde_p>B&&B_tilde_p<1+beta*B); % expenditure condition
    if (e_cond==1)&&(s_cond==1)
        [x,fval,flag]=gamma(B,B_tilde,B_tilde_p,params);
        R_out=1/beta*x;
    elseif (e_cond==1)&&(s_cond==0)
        [x,fval,flag]=eta(B,B_tilde,params);
        R_out=1/beta*x;
    elseif (e_cond==0)&&(s_cond==1)
        [x,fval,flag]=rho(B,B_tilde_p,params);
        R_out=1/beta*x;
    else 
        R_out=1/beta;
        fval=0; flag=1;
    end
end

function [Rg_out,fval,flag]=Rg(B,B_tilde,B_tilde_p,params)
    beta=params.beta;
    [Gamma,fval,flag]=gammaG(B,B_tilde,B_tilde_p,params);
    Rg_out=1/beta*Gamma;
end
