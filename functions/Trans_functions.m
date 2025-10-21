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

function []=plot_comp(index_pre,index_crunch,index_rec,index_post,var_s_ce,var_s_rp)
    % Colors: move to the top and pass as structure
    % Continuous time-type of plot
    color_t1=[0.7 0.7 0.7];
    color_t0=[0.0 0.0 0.6];
    linethick=3;

    % Plot pre-shock dates
    plot(index_pre,var_s_ce(index_pre),'LineWidth',linethick,'LineStyle','-','Color',color_t0); 
    scatter(index_pre(end),var_s_ce(index_pre(end)),'MarkerFaceColor','w','MarkerEdgeColor',color_t0); 

    % Plot fixed crunch phase
    plot(index_crunch-1,var_s_ce(index_crunch),'LineWidth',linethick,'LineStyle','-','Color',color_t0);
    plot(index_crunch-1,var_s_rp(index_crunch),'LineWidth',linethick,'LineStyle','--','Color',color_t1);
    scatter(index_crunch(1)-1,var_s_ce(index_crunch(1)),'filled','MarkerFaceColor',color_t0,'MarkerEdgeColor',color_t0);
    scatter(index_crunch(end)-1,var_s_ce(index_crunch(end)),'filled','MarkerFaceColor','w','MarkerEdgeColor',color_t0); 
    scatter(index_crunch(1)-1,var_s_rp(index_crunch(1)),'filled','MarkerFaceColor',color_t1,'MarkerEdgeColor',color_t1);
    scatter(index_crunch(end)-1,var_s_rp(index_crunch(end)),'filled','MarkerFaceColor','w','MarkerEdgeColor',color_t1); 


    % Plot recovery phase
    plot(index_rec-2,var_s_ce(index_rec),'LineWidth',linethick,'LineStyle','-','Color',color_t0);
    plot(index_rec-2,var_s_rp(index_rec),'LineWidth',linethick,'LineStyle','--','Color',color_t1);
    scatter(index_rec(1)-2,var_s_ce(index_rec(1)),'MarkerFaceColor',color_t0,'MarkerEdgeColor',color_t0);
    scatter(index_rec(end)-2,var_s_ce(index_rec(end)),'filled','MarkerFaceColor','w','MarkerEdgeColor',color_t0);
    scatter(index_rec(1)-2,var_s_rp(index_rec(1)),'MarkerFaceColor',color_t0,'MarkerEdgeColor',color_t0);
    scatter(index_rec(end)-2,var_s_rp(index_rec(end)),'filled','MarkerFaceColor','w','MarkerEdgeColor',color_t0);
    
    % Plot no shock phase
    plot(index_post-3,var_s_ce(index_post),'LineWidth',linethick,'LineStyle','-','Color',color_t0); 
    plot(index_post-3,var_s_rp(index_post),'LineWidth',linethick,'LineStyle','--','Color',color_t1);
    scatter(index_post(1)-3,var_s_ce(index_post(1)),'MarkerFaceColor','w','MarkerEdgeColor',color_t0); 
    scatter(index_post(end)-3,var_s_ce(index_post(end)),'>','filled','MarkerFaceColor',color_t0);
    scatter(index_post(1)-3,var_s_rp(index_post(1)),'MarkerFaceColor','w','MarkerEdgeColor',color_t0); 
    scatter(index_post(end)-3,var_s_rp(index_post(end)),'>','filled','MarkerFaceColor',color_t0);

    ybound=ylim;
    line([index_pre(end) index_pre(end)],ybound,'LineStyle',':','Color',[0.3 0.3 0.3]); 
    line([index_crunch(end)-1 index_crunch(end)-1],ybound,'LineStyle',':','Color',[0.3 0.3 0.3]); 
    line([index_rec(end)-2 index_rec(end)-2],ybound,'LineStyle',':','Color',[0.3 0.3 0.3]); 
    line([index_post(end)-3 index_post(end)-3],ybound,'LineStyle',':','Color',[0.3 0.3 0.3]);
    xlabel('Time','interpreter','latex');
end

function []=plot_comp2(index_pre,index_crunch,index_rec,index_post,var_a,var_b)
    % Color structures
    color_t1=[0.7 0.7 0.7];
    color_t0=[0.0 0.0 0.6];
    linethick=3;
    index_trans=[index_crunch index_rec index_post];

    % Plot pre-shock dates
    plot(index_pre,var_a(index_pre),'LineWidth',linethick,'LineStyle','-','Color',color_t0); 
    scatter(index_pre(end),var_a(index_pre(end)),'MarkerFaceColor','w','MarkerEdgeColor',color_t0); 

    % Plot fixed crunch phase
    plot(index_trans-1,var_a(index_trans),'LineWidth',linethick,'LineStyle','-','Color',color_t0);
    plot(index_trans-1,var_b(index_trans),'LineWidth',linethick,'LineStyle','--','Color',color_t1);
    scatter(index_trans(1)-1,var_a(index_trans(1)),'filled','MarkerFaceColor',color_t0,'MarkerEdgeColor',color_t0);
    scatter(index_trans(1)-1,var_b(index_trans(1)),'filled','MarkerFaceColor',color_t1,'MarkerEdgeColor',color_t1);

    % Plot no shock phase
    plot(index_trans-1,var_a(index_trans),'LineWidth',linethick,'LineStyle','-','Color',color_t0); 
    plot(index_trans-1,var_b(index_trans),'LineWidth',linethick,'LineStyle','--','Color',color_t1);
    scatter(index_trans(end)-1,var_a(index_trans(end)),'>','filled','MarkerFaceColor',color_t0);
    scatter(index_trans(end)-1,var_b(index_trans(end)),'>','filled','MarkerFaceColor',color_t1);

    ybound=ylim;
    line([index_pre(end) index_pre(end)],ybound,'LineStyle',':','Color',[0.3 0.3 0.3]); 
    line([index_crunch(end)-1 index_crunch(end)-1],ybound,'LineStyle',':','Color',[0.3 0.3 0.3]); 
    line([index_rec(end)-1 index_rec(end)-1],ybound,'LineStyle',':','Color',[0.3 0.3 0.3]); 
    xlabel('Time','interpreter','latex');
end