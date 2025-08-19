% Aggregate Euler Equations - Dynamics
figure('Name','Aggregate Euler')
range_test=[max(Bstar(B_tilde)-1,0)+0.9 Bast(B_tilde_p)+0.1];
fplot(@(x) Euler_t(x,B_tilde),range_test,'LineWidth',2,'Color',color_t0); hold on; grid on; % axis tight; 
fplot(@(x) Euler_t_p(x,B_tilde_p),range_test,'LineWidth',2,'LineStyle','--','Color',color_t1); title('Average Price'); hold on;
grid on; yvals=ylim; line([Bstar(B_tilde) Bstar(B_tilde)],[yvals(1) yvals(2)],'Color','k','LineWidth',1,'LineStyle','--');
grid on; yvals=ylim; line([B_ast_aux B_ast_aux],[yvals(1) yvals(2)],'Color','k','LineWidth',1,'LineStyle','--');
grid on; yvals=ylim; line([B_tilde B_tilde],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.5],'LineWidth',1,'LineStyle','-.'); grid minor;
grid on; yvals=ylim; line([Bhat(B_tilde,B_tilde) Bhat(B_tilde,B_tilde)],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.5],'LineWidth',1,'LineStyle',':');
title('Equilibrium Dynamics','interpreter', 'latex'); 
legend('$\mathcal{E}_{t}$','$\mathcal{E}_{t+1}$','Box','Off','interpreter', 'latex','Location','northwest','AutoUpdate','off')
set(gca,'XTick',Bstar(B_tilde),'XTickLabel',[]);
set(gca,'YTick',[yvals(1) yvals(2)],'YTickLabel',[]);
yTicks = get(gca,'ytick');
xTicks = get(gca, 'xtick');
ax = axis;    %Get left most x-position
text(ax(1)-0.01,yTicks(2),'$\mathcal{E}_t$', 'HorizontalAlignment','Right','interpreter', 'latex');
text(ax(1)-0.07,yTicks(2),'$\mathcal{E}_{t+1},$', 'HorizontalAlignment','Right','interpreter', 'latex');
text(ax(2)+0.2,yTicks(1)-0.1,'$B_{t},B_{t+1}$', 'HorizontalAlignment','Right','interpreter', 'latex');
text(xTicks(1),yTicks(1)-0.1,'$B^{\star}$', 'HorizontalAlignment','Center','VerticalAlignment','Middle','interpreter', 'latex');
text(Bast(B_tilde_p),yTicks(1)-0.1,'$B^{\ast}(\tilde{B})$', 'HorizontalAlignment','Center','VerticalAlignment','Middle','interpreter', 'latex');
text(B_tilde,yTicks(1)-0.1,'$\tilde{B}$', 'HorizontalAlignment','Center','VerticalAlignment','Middle','interpreter', 'latex'); hold on;
% add arrows of transition dynamics
drawArrow = @(x,y) quiver(x(1),y(1),x(2)-x(1),y(2)-y(1),0 ,'Color',[0.6290, 0.6940, 0.4250],'MaxHeadSize',0.45,'LineStyle','-','LineWidth',1,'Color',color_arrow) ;
B0_plot   = mean([B_tilde_p B_ast_aux]);
sim_plot  = sim_mod(B_tilde*ones(1,N_t),params,funcs,B0_plot);
B_p_t = sim_plot.B_t(2:end);
B_t   = sim_plot.B_t(1:end-1);
for ii = 1:N_t-1
    x=[B_t(ii) B_p_t(ii)];
    y=[Euler_t(B_t(ii),B_tilde) Euler_t_p(B_p_t(ii),B_tilde)];
    drawArrow(x,y); drawnow; pause(0.1)
    x=[B_p_t(ii) B_p_t(ii)];
    y=[Euler_t_p(B_p_t(ii),B_tilde) Euler_t(B_p_t(ii),B_tilde)];
    drawArrow(x,y); drawnow; pause(0.1)
end
if printit==1
    orient landscape;
    ax = gca;
    exportgraphics(ax,'F_transition.pdf','BackgroundColor','none');
end

%% Plot with Hysteresis Region
figure('Name','Aggregate Euler (Hysteresis)')
range_test=[max(Bstar(B_tilde)-1,0)+0.9 Bast(B_tilde_p)+1];
fplot(@(x) Euler_t(x,B_tilde),range_test,'LineWidth',2,'Color',color_t0); hold on; grid on; % axis tight; 
fplot(@(x) Euler_t_p(x,B_tilde_p),range_test,'LineWidth',2,'LineStyle','--','Color',color_t1); title('Average Price'); hold on;
grid on; yvals=ylim; line([Bstar(B_tilde) Bstar(B_tilde)],[yvals(1) yvals(2)],'Color','k','LineWidth',1,'LineStyle','--');
grid on; xvals=xlim; line([Bast(B_tilde_p) Bast(B_tilde_p)],[yvals(1) yvals(2)],'Color','k','LineWidth',1,'LineStyle','--');
grid on; yvals=ylim; line([B_tilde B_tilde],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.5],'LineWidth',1,'LineStyle','-.'); grid minor;
grid on; yvals=ylim; line([Bhat(B_tilde,B_tilde) Bhat(B_tilde,B_tilde)],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.5],'LineWidth',1,'LineStyle',':');
title('Equilibrium Dynamics','interpreter', 'latex'); 
legend('$\mathcal{E}_{t}$','$\mathcal{E}_{t+1}$','Box','Off','interpreter', 'latex','Location','northwest','AutoUpdate','off')
set(gca,'XTick',Bstar(B_tilde),'XTickLabel',[]);
set(gca,'YTick',[yvals(1) yvals(2)],'YTickLabel',[]);
yTicks = get(gca,'ytick');
xTicks = get(gca, 'xtick');
ax = axis; %Get left most x-position
text(ax(1)-0.01,yTicks(2),'$\mathcal{E}_t$', 'HorizontalAlignment','Right','interpreter', 'latex');
text(ax(1)-0.15,yTicks(2),'$\mathcal{E}_{t+1}$', 'HorizontalAlignment','Right','interpreter', 'latex');
text(ax(2),yTicks(1)-0.1,'$B_{t},B_{t+1}$', 'HorizontalAlignment','Right','interpreter', 'latex');
text(xTicks(1),yTicks(1)-0.1,'$B^{\star}$', 'HorizontalAlignment','Center','VerticalAlignment','Middle','interpreter', 'latex');
text(Bast(B_tilde_p),yTicks(1)-0.1,'$B^{\ast}(\tilde{B})$', 'HorizontalAlignment','Center','VerticalAlignment','Middle','interpreter', 'latex');
text(B_tilde,yTicks(1)-0.1,'$\tilde{B}$', 'HorizontalAlignment','Center','VerticalAlignment','Middle','interpreter', 'latex'); hold on;
x1 = [(Bast(B_tilde_p)+xvals(2))/2 (Bast(B_tilde_p)+xvals(2))/2 xvals(2) xvals(2)];
y1 = [yvals(1) yvals(2) yvals(2) yvals(1)];
fill(x1,y1,[0.5 0.5 0.5],'FaceAlpha',0.3);
if printit==1
    orient landscape;
     ax = gca;
    exportgraphics(ax,'F_Japan.pdf','BackgroundColor','none');
end

%% Outside of Steady States
figure('Name','Aggregate Euler')
range_test=[max(Bstar(B_tilde)-1,0) Bast(B_tilde_p)+1];
fplot(@(x) Euler_t(x,B_tilde),range_test,'LineWidth',2); title('RHS - LHS (extreme scenario)'); hold on; grid on; axis tight; 
fplot(@(x) Euler_t_p(x,B_tilde_p),range_test,'LineWidth',2,'LineStyle','--'); title('Average Price'); hold on;
grid on; yvals=ylim; line([Bstar(B_tilde) Bstar(B_tilde)],[yvals(1) yvals(2)],'Color','k','LineWidth',2,'LineStyle','--');
grid on; yvals=ylim; line([Bast(B_tilde_p) Bast(B_tilde_p)],[yvals(1) yvals(2)],'Color','k','LineWidth',2,'LineStyle','--');
grid on; yvals=ylim; line([B_tilde B_tilde],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.5],'LineWidth',2,'LineStyle','-.'); grid minor;
grid on; yvals=ylim; line([Bhat(B_tilde,B_tilde) Bhat(B_tilde,B_tilde)],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.5],'LineWidth',2,'LineStyle','-.');
legend('t Euler','t+1 Euler')

figure('Name','Equilibrium Plot')
fplot(@(x) B_p_res(B0,x,B_tilde,B_tilde_p),range_test,'LineWidth',2); title('Euler RHS residual Euler Equation');
grid on; yvals=ylim; line([Bstar(B_tilde) Bstar(B_tilde)],[yvals(1) yvals(2)],'Color','k','LineWidth',2,'LineStyle','--');
grid on; yvals=ylim; line([Bast(B_tilde) Bast(B_tilde)],[yvals(1) yvals(2)],'Color','k','LineWidth',2,'LineStyle','--');
grid on; yvals=ylim; line([B_tilde B_tilde],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.5],'LineWidth',2,'LineStyle','-.');
grid on; yvals=ylim; line([Bhat(B_tilde,B_tilde) Bhat(B_tilde,B_tilde)],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.5],'LineWidth',2,'LineStyle','-.');
grid on; yvals=ylim; line([Bhat(B_tilde,B_tilde) Bhat(B_tilde,B_tilde)],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.5],'LineWidth',2,'LineStyle','-.');
hold on;
plot(range_test,0*range_test,'LineWidth',2);
grid on; grid minor; axis tight;

