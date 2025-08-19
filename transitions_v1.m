%% Payment Crises: Code to Compute Equilibrium Solution 
% (c) Saki Bigio
clear; close all;

% Code Specs
printit=0; % Save plots into PDF
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
N_t=8;

% Baseline Version plots
if plotit==1
    Plots_AggregateEuler;
end

% Extreme version of plots
if plotit==2
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
plot(delta_vec,B_ast_vec); grid on;
xlabel('$\delta$','FontSize', 16,'Interpreter','latex');
ylabel('$B^{\ast}$','FontSize', 16,'Interpreter','latex');
load("ContinuousTimeSetting/bh_cutoffs_ct.mat")
hold on;
plot(delta_vec_ct,-b_h_delta_ct); grid on;
axis tight;

% Rest parameters
Trans_params; % reset parameters
Trans_function;

%% Steady-State in 3-d
% 3-d plots
if plotit==2
    Plots_SteadyState_3d;
end

%% Main Transition Simulation
B0_base     = 1.00              ; % initial value
B0          = B0_base           ; % Initial value for simuluations
B_tilde_rho = 0.6              ; % mean reversal rate after recovery
B_tilde_ss  = 0.05*params.B_bar ; % steady-state smooth borrowing limit
B_tilde_init= 0.1*B_tilde_ss   ; 
B_tilde_pre = B_tilde_ss        ;
B_star      = Bstar(B_tilde_ss) ;

% Example 1: Violent Crunch
T_pre      = 1  ; index_pre=1:T_pre;
T_crunch   = 20 ; index_crunch=T_pre+1:T_pre+T_crunch;
T_rec      = 100; index_rec=T_pre+T_crunch+1:T_pre+T_crunch+T_rec+1; 
T_post     = 10 ; index_post = T_pre + T_crunch + T_rec + 1 : T_pre + T_crunch + T_rec + T_post;

% Initial Values
B_tilde_t   = [ones(1,T_pre)*B_tilde_ss zeros(1,T_crunch+T_rec)+B_tilde_init ones(1,T_post)*B_tilde_ss]';
B_tilde_v   = B_tilde_t;

% Example #1: Violent Transitions
% T = length(B_tilde_t);
sim_violent=sim_mod(B_tilde_t,params,funcs,B0);

% Example #2: Smooth Transition  
% B_tilde_t   = [B_tilde_ss*ones(1,T_pre) zeros(1,T_crunch) B_tilde_ss-B_tilde_rho.^(1:T_rec)*B_tilde_ss ones(1,T_post)*B_tilde_ss]';
B_tilde_t   = [ones(1,T_pre)*B_tilde_ss zeros(1,T_crunch)+B_tilde_pre B_tilde_ss-B_tilde_rho.^(0:T_rec)*(B_tilde_ss-B_tilde_init) ones(1,T_post)*B_tilde_ss]';
B_tilde_s   = B_tilde_t;
T= length(B_tilde_t);
sim_smooth=sim_mod(B_tilde_s,params,funcs,B0);

% Plot Specs
varlist={'B_tilde_t','B_t','R_t','C_s_t','S_w_t','mu_t','q_t','X_w_t','Y_t'};
vartaglist={'$\tilde{B}_t$','$B_t$','$R_t$','$C^s_t$','$S^w_t$','$\mu_t$','$q_t$','$X^w_t$','$Y_t$'};
printlist={'fig_tildeB','fig_B','fig_R','fig_C','fig_S','fig_mu','fig_q','fig_X','fig_Y'};

% Plot of Smooth shock
if plotit==1
    for vv=1:numel(varlist)
      %  eval(['var_v=sim_violent.' varlist{vv} ';']);
        var_s=sim_smooth.(varlist{vv});
        vartag=vartaglist{vv};
        figure;
      %  plot_comp2(index_pre,index_crunch,index_rec,index_post,var_s,var_v_rp);
        title(vartag,'Interpreter','latex','FontSize',fontsize); hold on; grid on; 
        plot(index_pre,var_s(index_pre),'LineWidth',linethick,'LineStyle','-','Color',color_t0); hold on; 
        scatter(index_pre(end),var_s(index_pre(end)),'MarkerFaceColor','w','MarkerEdgeColor',color_t0); 
        plot(index_crunch-1,var_s(index_crunch),'LineWidth',linethick,'LineStyle','-','Color',color_t0); 
        scatter(index_crunch(end)-1,var_s(index_crunch(end)),'MarkerFaceColor','w','MarkerEdgeColor',color_t0); 
        scatter(index_crunch(end)-1,var_s(index_crunch(end)+1),'filled','MarkerFaceColor',color_t0,'MarkerEdgeColor',color_t0); 
        plot(index_rec-2,var_s(index_rec),'LineWidth',linethick,'LineStyle','-','Color',color_t0); 
        scatter(index_rec(end)-2,var_s(index_rec(end)),'MarkerFaceColor','w','MarkerEdgeColor',color_t0); 
        scatter(index_rec(end)-2,var_s(index_rec(end)),'MarkerFaceColor',color_t0,'MarkerEdgeColor',color_t0);
        plot(index_post-3,var_s(index_post),'LineWidth',linethick,'LineStyle','-','Color',color_t0); grid on; % axis tight;
        scatter(index_post(end)-3,var_s(index_post(end)),'MarkerFaceColor','w','MarkerEdgeColor',color_t0); 
        scatter(index_post(end)-3,var_s(index_post(end)),'MarkerFaceColor',color_t0,'MarkerEdgeColor',color_t0);
        scatter(index_post(end)-3,var_s(index_post(end)),'>','filled','MarkerFaceColor',color_t0); 
        ybound=ylim;
        line([index_pre(end) index_pre(end)],ybound,'LineStyle',':','Color',[0.3 0.3 0.3]); 
        line([index_crunch(end)-1 index_crunch(end)-1],ybound,'LineStyle',':','Color',[0.3 0.3 0.3]); 
        line([index_rec(end)-2 index_rec(end)-2],ybound,'LineStyle',':','Color',[0.3 0.3 0.3]); 
        line([index_post(end)-3 index_post(end)-3],ybound,'LineStyle',':','Color',[0.3 0.3 0.3]);
        xlabel('Time','interpreter','latex');
        if 1>1
            legend('Violent Transition','Smooth Transition','Interpreter','Latex','FontSize',fontsize,'Box','off','Location','NorthWest');
        end
        orient landscape;
        eval(['saveas(gcf,''' printlist{vv} ''',''pdf'');']);
    end
end

%% Resilience Plot
% run same plot with lower delta
params.delta=0.85;
Trans_function; % Important, update functions
B0          = B0_base;
sim_deltachange=sim_mod(B_tilde_s,params,funcs,B0);
for vv=1:numel(varlist)
    var_baseline=sim_smooth.(varlist{vv});
    var_higherd =sim_deltachange.(varlist{vv});
    vartag=vartaglist{vv};
    figure('Name','Resilience Plot');
    title(vartag,'Interpreter','latex','FontSize',fontsize); hold on; grid on; 
    plot_comp2(index_pre,index_crunch,index_rec,index_post,var_baseline,var_higherd);
    if vv==1
        legend('High Delta','Lower Delta','Interpreter','Latex','FontSize',fontsize,'Location','NorthWest','Box','off');
    end
    if printit==1
        orient landscape;
    end
end

%% Planner Solutions
B_ss_violent=sim_violent.B_t(end);
theta_violent=1-(1-params.beta)*B_ss_violent;
sim_p_violent=sim_planner(B_tilde_v,theta_violent,params,funcs,B0);

% Main Simulation
printlist={'fig2_tildeB','fig2_B','fig2_R','fig2_C','fig2_S','fig2_mu','fig2_q','fig2_X','fig2_Y'};
if plotit==1
    for vv=1:numel(varlist)
        var_v_ce=sim_violent.(varlist{vv});
        var_v_rp=sim_p_violent.(varlist{vv});
        vartag=vartaglist{vv};
        plot_comp2(index_pre,index_crunch,index_rec,index_post,var_v_ce,var_v_rp);
        title(vartag,'Interpreter','latex','FontSize',fontsize); hold on; grid on; 
        xlabel('Time','Interpreter','latex','FontSize',fontsize); title(vartag,'Interpreter','latex','FontSize',fontsize); hold on; 
        if vv==1
            legend('Violent Transition (CE)','Violent Transition (PE)','Interpreter','Latex','FontSize',fontsize,'Location','NorthWest','Box','off');
        end
        orient landscape;
        eval(['saveas(gcf,''' printlist{vv} ''',''pdf'');']);
    end
end

% Planner vs. Market - Smooth Transition
B_ss_smooth=sim_smooth.B_t(end);
theta_smooth=1-(1-params.beta)*B_ss_smooth;
sim_p_smooth=sim_planner(B_tilde_s,theta_smooth,params,funcs,B0);
printlist={'fig3_tildeB','fig3_B','fig3_R','fig3_C','fig3_S','fig3_mu','fig3_q','fig3_X','fig3_Y'};
for vv=1:numel(varlist)
    var_s_ce=sim_smooth.(varlist{vv});
    var_s_rp=sim_p_smooth.(varlist{vv});
    vartag=vartaglist{vv};
    figure('Name','Comp-Planner (smooth shock)');
    title(vartag,'Interpreter','latex','FontSize',fontsize); hold on; grid on; 
    plot_comp2(index_pre,index_crunch,index_rec,index_post,var_s_ce,var_s_rp);
    if vv==1
        legend('Competitive Equilbrium','Planner Solution','Interpreter','Latex','FontSize',fontsize,'Location','NorthWest','Box','off');
    end
    if printit==1
        orient landscape;
        eval(['saveas(gcf,''' printlist{vv} ''',''pdf'');']);
    end
end

%% Calling Functions
% Transition plot functions
Trans_Plot_functions;