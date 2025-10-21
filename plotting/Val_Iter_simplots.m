function []=Val_Iter_simplots(hist,params,B_exit_num)
% print options
printit=params.printit;

% Color Specs
color_t0=[0.0 0.0 0.6];
linethick=3;
fontsize=16;

% Color Ordering for Plots
newcolors = [0.00 0.00 1.0
             0.5 0.5 0.9
             0.25 0.75 0.9
             0.25 0.20 0.9
             0.7 0.7 0.7];
marker_size=120;

% Updating
[N_t,N_c]=size(hist.B_hist);
time=1:N_t;

% Plot Specs
varlist={'B_hist','C_hist','X_hist','S_hist'};
vartaglist={'$B_t$','$C_t$','$X_t$','$S_t$'};
printlist={'F_path_B','fig_path_C','fig_path_S','fig_path_X'};
for vv=1:numel(varlist)
    var_s=hist.(varlist{vv});
    vartag=vartaglist{vv};
    
    figure('Name',string(vartag));
   %  title(string(vartag),'Interpreter','latex','FontSize',fontsize); 
    hold on; grid on; 
    plot(time,var_s,'LineWidth',linethick); 
    linestyleorder("mixedstyles");
    colororder(newcolors);
    scatter(ones(1,N_c),var_s(1,1:end),marker_size,'MarkerFaceColor','w','MarkerEdgeColor',color_t0);
    scatter(time(end),var_s(end,1:end),marker_size,'>','filled','MarkerFaceColor','w','MarkerEdgeColor','k');
    xlim([time(1) time(end)]);
   % xlabel('Time (t)','interpreter','latex');
    if vv==1
        yvals=ylim; xvals=xlim;
        set(gca,'YTick',params.Bstar,'YTickLabel',[]);
   %     set(gca,'YTick',params.Bast,'YTickLabel',[]);
        set(gca,'YTick',params.B_tilde,'YTickLabel',[]);
        set(gca,'YTick',params.b_bar,'YTickLabel',[]);
        set(gca,'YTick',B_exit_num,'YTickLabel',[]);
        set(gca,'XTick',[xvals(1) xvals(2)],'XTickLabel',[]);
        yTicks = get(gca,'ytick');
        xTicks = get(gca,'xtick');
        text(xTicks(1)-0.2,params.Bstar,'$B^{\star}$', 'HorizontalAlignment','Center','VerticalAlignment','top','interpreter', 'latex','FontSize',14);
    %    text(xTicks(1)-1,params.Bast,'$B^{\ast}(\tilde{B})$', 'HorizontalAlignment','Center','VerticalAlignment','top','interpreter', 'latex','FontSize',14);
        text(xTicks(1)-0.2,params.B_tilde,'$\tilde{B}$', 'HorizontalAlignment','Center','VerticalAlignment','top','interpreter', 'latex','FontSize',14); hold on;
        text(xTicks(1)-0.2,params.b_bar,'$\bar{B}$', 'HorizontalAlignment','Center','VerticalAlignment','top','interpreter', 'latex','FontSize',14); hold on;
        text(xTicks(1)-0.2,B_exit_num,'$B^{h}$', 'HorizontalAlignment','Center','VerticalAlignment','top','interpreter', 'latex','FontSize',14); hold on;
        text(xvals(2)+0.2,yvals(1)-0.1,'Time (t)', 'HorizontalAlignment','Center','VerticalAlignment','top','interpreter', 'latex','FontSize',14); hold on;
        text(xvals(1)-0.2,yvals(2),'$B_t$', 'HorizontalAlignment','Center','VerticalAlignment','top','interpreter', 'latex','FontSize',14); hold on;

        line([time(1) time(end)],[B_exit_num B_exit_num],'Color','k','LineWidth',1,'LineStyle','--');
        line([time(1) time(end)],[params.B_tilde params.B_tilde],'Color','k','LineWidth',1,'LineStyle','--');
    %    line([time(1) time(end)],[params.Bast params.Bast],'Color','k','LineWidth',1,'LineStyle','--');
        line([time(1) time(end)],[params.Bstar params.Bstar],'Color','k','LineWidth',1,'LineStyle','--');

        if printit==1    
                orient landscape;
                ax = gca;
                exportgraphics(ax,[params.folder printlist{vv} '.pdf'],'BackgroundColor','none');
        end
    end
%    if printit==1
 %       orient landscape;
 %       eval(['saveas(gcf,''' printlist{vv} ''',''pdf'');']);
 %   end
end