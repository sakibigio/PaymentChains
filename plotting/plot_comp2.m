function []=plot_comp2(index_pre,index_crunch,index_rec,index_post,var_a,var_b)
    % Color structures
    color_t0=[0.0 0.0 0.6];
    color_t1=[0.7 0.7 0.7];
    linethick=3;
    index_trans=[index_pre index_crunch index_rec index_post];

    % Plot pre-shock dates
    plot(index_pre,var_a(index_pre),'LineWidth',linethick,'LineStyle','-','Color',color_t0); 
    plot(index_pre,var_b(index_pre),'LineWidth',linethick,'LineStyle','--','Color',color_t1); 
    scatter(index_pre(end)-1,var_a(index_pre(end)),'MarkerFaceColor','w','MarkerEdgeColor',color_t0); 

    % Plot fixed crunch phase
    plot(index_trans-1,var_a(index_trans),'LineWidth',linethick,'LineStyle','-','Color',color_t0);
    plot(index_trans-1,var_b(index_trans),'LineWidth',linethick,'LineStyle','--','Color',color_t1);
    scatter(index_rec(1)-1,var_a(index_rec(1)),40,'filled','MarkerFaceColor',color_t0,'MarkerEdgeColor',color_t0);
    scatter(index_rec(1)-1,var_b(index_rec(1)),40,'filled','MarkerFaceColor',color_t1,'MarkerEdgeColor',color_t1);

    % Plot no shock phase
%    plot(index_trans-1,var_a(index_trans),'LineWidth',linethick,'LineStyle','-','Color',color_t0); 
%    plot(index_trans-1,var_b(index_trans),'LineWidth',linethick,'LineStyle','--','Color',color_t1);
    scatter(index_trans(end)-1,var_a(index_trans(end)),'>','filled','MarkerFaceColor',color_t0);
    scatter(index_trans(end)-1,var_b(index_trans(end)),'>','filled','MarkerFaceColor',color_t1);

    ybound=ylim;
    % line([index_pre(end) index_pre(end)],ybound,'LineStyle',':','Color',[0.3 0.3 0.3]); 
    line([index_crunch(end)-1 index_crunch(end)-1],ybound,'LineStyle',':','Color',[0.3 0.3 0.3],'LineWidth',2); 
    % line([index_rec(end)-1 index_rec(end)-1],ybound,'LineStyle',':','Color',[0.3 0.3 0.3]); 
    xlabel('Time','interpreter','latex');
end