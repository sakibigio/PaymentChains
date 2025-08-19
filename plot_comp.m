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