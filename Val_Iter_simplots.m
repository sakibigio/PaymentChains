function []=Val_Iter_simplots(hist)
% print options
printit=0;

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
printlist={'fig_path_B','fig_path_C','fig_path_S','fig_path_X'};
for vv=1:numel(varlist)
    var_s=hist.(varlist{vv});
    vartag=vartaglist{vv};
    
    figure;
    title(string(vartag),'Interpreter','latex','FontSize',fontsize); hold on; grid on; 
    plot(time,var_s,'LineWidth',linethick); 
    linestyleorder("mixedstyles");
    colororder(newcolors);
    scatter(ones(1,N_c),var_s(1,1:end),marker_size,'MarkerFaceColor','w','MarkerEdgeColor',color_t0);
    scatter(time(end),var_s(end,1:end),marker_size,'>','filled','MarkerFaceColor','w','MarkerEdgeColor','k');
    xlim([time(1) time(end)]);
    xlabel('Time (t)','interpreter','latex');
    if printit==1
        orient landscape;
        eval(['saveas(gcf,''' printlist{vv} ''',''pdf'');']);
    end
end