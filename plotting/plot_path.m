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