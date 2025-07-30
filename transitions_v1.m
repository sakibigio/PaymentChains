%% Payment Crises: Code to Compute Equilibrium Solution 
% (c) Saki Bigio
clear; close all;

% Code Specs
printit=1; % Save plots into PDF
plotit=0; % 1 if plot needed

%% Parameter Set
delta=0.9; % Value of discount in production
beta =0.95; % Discount factor
B_bar=1/(1-beta); 
fontsize=16;

% Plot Specs
color_t1=[0.7 0.7 0.7];
color_t0=[0.0 0.0 0.6];
color_arrow=[0 0 0];

%% Equilibrium Functions
% CalA_ex=@(mu) delta/(1-delta)*(1-mu)./mu.*log((1-mu.*delta)./(1-mu));
if plotit==1
    CalA_ex=@(mu,delta) delta/(1-delta)*(1-mu)./mu.*log((1-mu.*delta)./(1-mu));
    figure
    fplot(@(x) CalA_ex(x,0.9),[0.01 0.99],'LineWidth',2,'Color',color_t0); hold on;
    fplot(@(x) CalA_ex(x,0.5),[0.01 0.99],'LineWidth',2,'Color',color_t1); xlabel('$\mu$','Interpreter','latex','FontSize',fontsize);
    title('$\mathcal{A}(\mu)$','interpreter','latex','FontSize',fontsize); grid on; axis tight; hold on;
    legend('High $\delta$','Low $\delta$','Interpreter','Latex')
    if printit==1
        orient landscape;
        ax = gca;
        exportgraphics(ax,'TFP_example.pdf','BackgroundColor','none');
        % saveas(gcf, 'TFP_example', 'pdf')
    end

    figure;
    chain=0:1:8;
    mu1=0.6; mu2=0.7;
    b=bar(chain,[(1-mu1)*mu1.^(chain);(1-mu2)*mu1.^(chain)],'BarWidth',1.2); hold on;
    b(1).FaceColor=color_t0;
    b(2).FaceColor=color_t1; xlabel('Chain Length','Interpreter','latex');
    title('Probability Distribution','Interpreter','latex'); grid on; axis tight; 
    legend('Low $\mu$','High $\mu$','Interpreter','Latex');
     if printit==1
        orient landscape;
        % saveas(gcf, 'Dist_example', 'pdf')
        ax = gca;
        exportgraphics(ax,'Dist_example.pdf','BackgroundColor','none');

    end
end



% q_mu=@(mu) CalA(mu);
if plotit==1
    figure
    fplot(q_mu,[0.01 0.99]); xlabel('$\mu$','Interpreter','latex');
    title('$q(\mu)$','Interpreter','latex'); grid on; axis tight;
end

% Main Functional forms
C_s   = @(B) (1-beta).*B;
S_w   = @(B,B_tilde) min(max(B_tilde-B,0),1-(1-beta)*B);
mu_Z  = @(B,B_tilde) 1-C_s(B)-S_w(B,B_tilde);
q_Z   = @(B,B_tilde) q_mu(mu_Z(B,B_tilde),delta);
X_w   = @(B,B_tilde) mu_Z(B,B_tilde)./q_Z(B,B_tilde);
Y_Z   = @(B,B_tilde) 1-mu_Z(B,B_tilde)+mu_Z(B,B_tilde)./q_mu(mu_Z(B,B_tilde),delta);

% Derived Functions
% Average price
Q_Z   = @(B,B_tilde) (1-(1-beta)*B)./(X_w(B,B_tilde)+S_w(B,B_tilde));

% Marginal Price of Savings
q_B=@(B,B_tilde) (1+(q_Z(B,B_tilde)-1).*(B_tilde<B));

% Marginal Price of eExpenditures
q_E=@(B,B_tilde) (1+(q_Z(B,B_tilde)-1).*(B_tilde<1+beta*B));

% Marginal Inflation
Pi =@(B,B_p,B_tilde,B_tilde_p) q_B(B_p,B_tilde_p)./q_E(B,B_tilde);

% Ratio of Average to Marginal Price (expenditure and savings)
Qq_E=@(B,B_tilde) Q_Z(B,B_tilde)./q_E(B,B_tilde);
Qq_B=@(B,B_tilde) Q_Z(B,B_tilde)./q_B(B,B_tilde);

% Euler_rhs=@(B,B_p,B_tilde,B_tilde_p)  B .* (1-beta*B_p)./(1-beta*B).*Q_Z(B,B_tilde)./Q_Z(B_p,B_tilde_p).*Pi(B,B_p,B_tilde,B_tilde_p);
% Euler_rhs=@(B,B_p,B_tilde,B_tilde_p)  B./(1-(1-beta)*B).*Q_Z(B,B_tilde).*(X_w(B_p,B_tilde_p)+S_w(B_p,B_tilde_p)).*Pi(B,B_p,B_tilde,B_tilde_p);
Euler_t_p = @(B_p,B_tilde_p) 1./(1/B_p-(1-beta)).*Qq_B(B_p,B_tilde_p);
Euler_t   = @(B,B_tilde) 1./(1/B-(1-beta)).*Qq_E(B,B_tilde);

B_p_res=@(B,B_p,B_tilde,B_tilde_p) Euler_t_p(B_p,B_tilde_p)  - Euler_t(B,B_tilde);

%B_p_res=@(B,B_p,B_tilde,B_tilde_p) (1/B-(1-beta))./Qq_E(B,B_tilde)-(1/B_tilde-(1-beta))./Qq_B(B_tilde,B_tilde_p);

% Critical Points
Bstar = @(B_tilde) 1/beta*(B_tilde-1);
Bast  = @(B_tilde_p) B_tilde_p/(q_Z(B_tilde_p,B_tilde_p).^(-1)-(q_Z(B_tilde_p,B_tilde_p).^(-1)-1)*(1-beta)*B_tilde_p);
Bhat  = @(B_tilde,B_tilde_p) fsolve(@(x) Euler_t(x,B_tilde)-1./(1/Bstar(B_tilde)-(1-beta)),Bstar(B_tilde)+0.2); 

%% Key solutions for planner
P    = @(B,B_tilde,theta) (1-theta)*log((1-beta)*B)+theta*log(S_w(B,B_tilde)+X_w(B,B_tilde));
B_Ram= @(B_tilde,theta) fminsearch(@(x) -P(x,B_tilde,theta),0.01);

%% Checks: all satisfied
% make sure C+S+X adds to Y
y_check=@(B,B_tilde) Y_Z(B,B_tilde)-(C_s(B)+S_w(B,B_tilde)+X_w(B,B_tilde));

%% Saving Functions
% make sure C+S+q^(-1)X adds to 1
income_check=@(B,B_tilde) 1-(C_s(B)+S_w(B,B_tilde)+q_Z(B,B_tilde).*X_w(B,B_tilde));
params.y_check=y_check;
params.income_check=income_check;
params.C_s=C_s;
params.S_w=S_w;
params.mu_Z=mu_Z;
params.q_Z=q_Z;
params.X_w=X_w;
params.Y_Z=Y_Z;
params.q_Z=q_Z;
params.beta=beta; 
params.delta=delta;
params.Q_Z=Q_Z;
params.Qq_E=Qq_E;
params.Qq_B=Qq_B;
params.Bstar=Bstar;
params.Bast=Bast;
params.B_p_res=B_p_res;

params.B_Ram=B_Ram;
params.P=P;

%% Test Plots
B0=1.5;
range_test=[0.9 3]; 

B_tilde=2; B_tilde_p=2;
figure
fplot(@(x) q_E(x,B_tilde),range_test,'LineWidth',2); title('Borrowing Marginal Price'); axis tight;
grid on; yvals=ylim; line([B_tilde B_tilde],[yvals(1) yvals(2)],'Color','k','LineWidth',2);

figure
fplot(@(x) q_B(x,B_tilde),range_test,'LineWidth',2); title('Borrowing Marginal Price');
grid on; yvals=ylim; line([B_tilde B_tilde],[yvals(1) yvals(2)],'Color','k','LineWidth',2);

figure
fplot(@(x) Q_Z(x,B_tilde_p),range_test,'LineWidth',2); title('Average Price');
grid on; yvals=ylim; line([B_tilde B_tilde],[yvals(1) yvals(2)],'Color','k','LineWidth',2);

figure
fplot(@(x) Q_Z(x,B_tilde_p)./Q_Z(B0,B_tilde),range_test,'LineWidth',2); title('Ratio of Average Prices Ratio');
grid on; yvals=ylim; line([B_tilde B_tilde],[yvals(1) yvals(2)],'Color','k','LineWidth',2);

figure
fplot(@(x) Pi(2,x,B_tilde,B_tilde_p),range_test,'LineWidth',2); title('Marginal Inflation');
grid on; yvals=ylim; line([B_tilde B_tilde],[yvals(1) yvals(2)],'Color','k','LineWidth',2);

% figure
% fplot(@(x) -Euler_rhs(B0,x,B_tilde,B_tilde_p),range_test,'LineWidth',2); title(' Euler RHS (must be dec)');
% hold on;
% plot(range_test,range_test,'LineWidth',2,'Color','k','LineStyle','--');
% grid on; axis tight;


%% Transitional Dynamics
figure
B_tilde=2; B_tilde_p=B_tilde;
range_test=[max(Bstar(B_tilde)-1,0)+0.9 Bast(B_tilde_p)+0.1];
fplot(@(x) Euler_t(x,B_tilde),range_test,'LineWidth',2,'Color',color_t0); hold on; grid on; % axis tight; 
fplot(@(x) Euler_t_p(x,B_tilde_p),range_test,'LineWidth',2,'LineStyle','--','Color',color_t1); title('Average Price'); hold on;
grid on; yvals=ylim; line([Bstar(B_tilde) Bstar(B_tilde)],[yvals(1) yvals(2)],'Color','k','LineWidth',1,'LineStyle','--');
grid on; yvals=ylim; line([Bast(B_tilde_p) Bast(B_tilde_p)],[yvals(1) yvals(2)],'Color','k','LineWidth',1,'LineStyle','--');
grid on; yvals=ylim; line([B_tilde B_tilde],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.5],'LineWidth',1,'LineStyle','-.'); grid minor;
grid on; yvals=ylim; line([Bhat(B_tilde,B_tilde) Bhat(B_tilde,B_tilde)],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.5],'LineWidth',1,'LineStyle',':');
title('Equilibrium Dynamics','interpreter', 'latex'); 
legend('$\mathcal{E}_{t}$','$\mathcal{E}_{t+1}$','Box','Off','interpreter', 'latex','Location','northwest','AutoUpdate','off')
% grid on; yvals=ylim; xvals=xlim; line([xvals(1) xvals(end)],[1/beta 1/beta],'Color','r','LineWidth',2,'LineStyle',':');

B_star_aux=Bstar(B_tilde);
% axis([0 1.9 0 1.5]);
% set the ticks and labels
set(gca,'XTick',[B_star_aux],'XTickLabel',[]);
set(gca,'YTick',[yvals(1) yvals(2)],'YTickLabel',[]);
yTicks = get(gca,'ytick');
xTicks = get(gca, 'xtick');

ax = axis; %Get left most x-position
text(ax(1)-0.01,yTicks(2),'$\mathcal{E}_t$', 'HorizontalAlignment','Right','interpreter', 'latex');
text(ax(1)-0.07,yTicks(2),'$\mathcal{E}_{t+1},$', 'HorizontalAlignment','Right','interpreter', 'latex');
text(ax(2)+0.2,yTicks(1)-0.1,'$B_{t},B_{t+1}$', 'HorizontalAlignment','Right','interpreter', 'latex');
% text(ax(1)-0.01,1.45,'$y,i,\delta k$', 'HorizontalAlignment','Right','interpreter', 'latex');
text(xTicks(1),yTicks(1)-0.1,'$B^{\star}$', 'HorizontalAlignment','Center','VerticalAlignment','Middle','interpreter', 'latex');
text(Bast(B_tilde_p),yTicks(1)-0.1,'$B^{\ast}(\tilde{B})$', 'HorizontalAlignment','Center','VerticalAlignment','Middle','interpreter', 'latex');
text(B_tilde,yTicks(1)-0.1,'$\tilde{B}$', 'HorizontalAlignment','Center','VerticalAlignment','Middle','interpreter', 'latex'); hold on;
% text(1.85,-0.03,'$k$', 'HorizontalAlignment','Center','VerticalAlignment','Middle','interpreter', 'latex');

% Main Iterations
% add arrows of transition dynamics
N_t=8;
drawArrow = @(x,y) quiver(x(1),y(1),x(2)-x(1),y(2)-y(1),0 ,'Color',[0.6290, 0.6940, 0.4250],'MaxHeadSize',0.45,'LineStyle','-','LineWidth',1,'Color',color_arrow) ;
B0_plot   = mean([B_tilde_p Bast(B_tilde_p)]);
sim_plot  = sim_mod(B_tilde*ones(1,N_t),params,B0_plot);

B_p_t = sim_plot.B_t(2:end);
B_t   = sim_plot.B_t(1:end-1);
for ii = 1:N_t-1
    x=[B_t(ii) B_p_t(ii)];
    y=[Euler_t(B_t(ii),B_tilde) Euler_t_p(B_p_t(ii),B_tilde)];
    drawArrow(x,y); drawnow; pause(0.1)
    x=[B_p_t(ii) B_p_t(ii)];
    y=[Euler_t_p(B_p_t(ii),B_tilde) Euler_t(B_p_t(ii),B_tilde)];
    drawArrow(x,y); drawnow; pause(1)
end
if printit==1
    orient landscape;
    % saveas(gcf,'F_transition','pdf');
    ax = gca;
    exportgraphics(ax,'F_transition.pdf','BackgroundColor','none');
end

% added lines for purpose of presentation
% drawbrace([1 syss],[1 yss],10, 'Color', 'k');

% legends
% text(max(k)+0.05,max(y),'$y$','interpreter', 'latex');
% text(max(k)+0.05,max(sy),'$sy$','interpreter', 'latex');
% text(max(k)+0.05,max(de),'$\delta k$','interpreter', 'latex');
 % text(kss-0.27,(syss+yss)/2,'$c=(1-s)y$','interpreter', 'latex');

%% Japan Steady State
figure
B_tilde=2; B_tilde_p=B_tilde;
range_test=[max(Bstar(B_tilde)-1,0)+0.9 Bast(B_tilde_p)+1];
fplot(@(x) Euler_t(x,B_tilde),range_test,'LineWidth',2,'Color',color_t0); hold on; grid on; % axis tight; 
fplot(@(x) Euler_t_p(x,B_tilde_p),range_test,'LineWidth',2,'LineStyle','--','Color',color_t1); title('Average Price'); hold on;
grid on; yvals=ylim; line([Bstar(B_tilde) Bstar(B_tilde)],[yvals(1) yvals(2)],'Color','k','LineWidth',1,'LineStyle','--');
grid on; xvals=xlim; line([Bast(B_tilde_p) Bast(B_tilde_p)],[yvals(1) yvals(2)],'Color','k','LineWidth',1,'LineStyle','--');
grid on; yvals=ylim; line([B_tilde B_tilde],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.5],'LineWidth',1,'LineStyle','-.'); grid minor;
grid on; yvals=ylim; line([Bhat(B_tilde,B_tilde) Bhat(B_tilde,B_tilde)],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.5],'LineWidth',1,'LineStyle',':');
title('Equilibrium Dynamics','interpreter', 'latex'); 
legend('$\mathcal{E}_{t}$','$\mathcal{E}_{t+1}$','Box','Off','interpreter', 'latex','Location','northwest','AutoUpdate','off')
% grid on; yvals=ylim; xvals=xlim; line([xvals(1) xvals(end)],[1/beta 1/beta],'Color','r','LineWidth',2,'LineStyle',':');

B_star_aux=Bstar(B_tilde);
% axis([0 1.9 0 1.5]);
% set the ticks and labels
set(gca,'XTick',B_star_aux,'XTickLabel',[]);
set(gca,'YTick',[yvals(1) yvals(2)],'YTickLabel',[]);
yTicks = get(gca,'ytick');
xTicks = get(gca, 'xtick');

ax = axis; %Get left most x-position
text(ax(1)-0.01,yTicks(2),'$\mathcal{E}_t$', 'HorizontalAlignment','Right','interpreter', 'latex');
text(ax(1)-0.15,yTicks(2),'$\mathcal{E}_{t+1}$', 'HorizontalAlignment','Right','interpreter', 'latex');
text(ax(2),yTicks(1)-0.1,'$B_{t},B_{t+1}$', 'HorizontalAlignment','Right','interpreter', 'latex');
% text(ax(1)-0.01,1.45,'$y,i,\delta k$', 'HorizontalAlignment','Right','interpreter', 'latex');
text(xTicks(1),yTicks(1)-0.1,'$B^{\star}$', 'HorizontalAlignment','Center','VerticalAlignment','Middle','interpreter', 'latex');
text(Bast(B_tilde_p),yTicks(1)-0.1,'$B^{\ast}(\tilde{B})$', 'HorizontalAlignment','Center','VerticalAlignment','Middle','interpreter', 'latex');
text(B_tilde,yTicks(1)-0.1,'$\tilde{B}$', 'HorizontalAlignment','Center','VerticalAlignment','Middle','interpreter', 'latex'); hold on;
% text(1.85,-0.03,'$k$', 'HorizontalAlignment','Center','VerticalAlignment','Middle','interpreter', 'latex');
x1 = [(Bast(B_tilde_p)+xvals(2))/2 (Bast(B_tilde_p)+xvals(2))/2 xvals(2) xvals(2)];
y1 = [yvals(1) yvals(2) yvals(2) yvals(1)];
fill(x1,y1,[0.5 0.5 0.5],'FaceAlpha',0.5);

% Main Iterations
% add arrows of transition dynamics
% N_t=8;
% drawArrow = @(x,y) quiver(x(1),y(1),x(2)-x(1),y(2)-y(1),0 ,'Color',[0.6290, 0.6940, 0.4250],'MaxHeadSize',0.45,'LineStyle','-','LineWidth',1) ;
% B0_plot   = mean([B_tilde_p Bast(B_tilde_p)]);
% sim_plot  = sim_mod(B_tilde*ones(1,N_t),params,B0_plot);
% 
% B_p_t = sim_plot.B_t(2:end);
% B_t   = sim_plot.B_t(1:end-1);
% for ii = 1:N_t-1
%     x=[B_t(ii) B_p_t(ii)];
%     y=[Euler_t(B_t(ii),B_tilde) Euler_t_p(B_p_t(ii),B_tilde)];
%     drawArrow(x,y); drawnow; pause(1)
%     x=[B_p_t(ii) B_p_t(ii)];
%     y=[Euler_t_p(B_p_t(ii),B_tilde) Euler_t(B_p_t(ii),B_tilde)];
%     drawArrow(x,y); drawnow; pause(1)
% end
if printit==1
    orient landscape;
   %  saveas(gcf,'F_Japan','pdf');
     ax = gca;
    exportgraphics(ax,'F_Japan.pdf','BackgroundColor','none');
end

%% Outside of Steady States

figure
B_tilde=2; B_tilde_p=2;
range_test=[max(Bstar(B_tilde)-1,0) Bast(B_tilde_p)+1];
fplot(@(x) Euler_t(x,B_tilde),range_test,'LineWidth',2); title('RHS - LHS (extreme scenario)'); hold on; grid on; axis tight; 
fplot(@(x) Euler_t_p(x,B_tilde_p),range_test,'LineWidth',2,'LineStyle','--'); title('Average Price'); hold on;
grid on; yvals=ylim; line([Bstar(B_tilde) Bstar(B_tilde)],[yvals(1) yvals(2)],'Color','k','LineWidth',2,'LineStyle','--');
grid on; yvals=ylim; line([Bast(B_tilde_p) Bast(B_tilde_p)],[yvals(1) yvals(2)],'Color','k','LineWidth',2,'LineStyle','--');
grid on; yvals=ylim; line([B_tilde B_tilde],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.5],'LineWidth',2,'LineStyle','-.'); grid minor;
grid on; yvals=ylim; line([Bhat(B_tilde,B_tilde) Bhat(B_tilde,B_tilde)],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.5],'LineWidth',2,'LineStyle','-.');
legend('t Euler','t+1 Euler')


figure
fplot(@(x) B_p_res(B0,x,B_tilde,B_tilde_p),range_test,'LineWidth',2); title('Euler RHS residual Euler Equation');
grid on; yvals=ylim; line([Bstar(B_tilde) Bstar(B_tilde)],[yvals(1) yvals(2)],'Color','k','LineWidth',2,'LineStyle','--');
grid on; yvals=ylim; line([Bast(B_tilde) Bast(B_tilde)],[yvals(1) yvals(2)],'Color','k','LineWidth',2,'LineStyle','--');
grid on; yvals=ylim; line([B_tilde B_tilde],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.5],'LineWidth',2,'LineStyle','-.');
grid on; yvals=ylim; line([Bhat(B_tilde,B_tilde) Bhat(B_tilde,B_tilde)],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.5],'LineWidth',2,'LineStyle','-.');
grid on; yvals=ylim; line([Bhat(B_tilde,B_tilde) Bhat(B_tilde,B_tilde)],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.5],'LineWidth',2,'LineStyle','-.');

% grid on; yvals=ylim; xvals=xlim; line([xvals(1) xvals(end)],[1/beta 1/beta],'Color','r','LineWidth',2,'LineStyle',':');
hold on;

plot(range_test,0*range_test,'LineWidth',2);
grid on; grid minor; axis tight;

%% Steady-State Plots
B_tilde= B_bar/5;
B_vec  = linspace(Bstar(B_tilde)-1,Bast(B_tilde)+2,2000)  ;
R_vec  = NaN*B_vec; B_p_vec  = NaN*B_vec; val_vec  = NaN*B_vec; flag_vec  = NaN*B_vec;
ii=0;
for B0=B_vec
    ii=ii+1;
    [B_p_vec(ii),val_vec(ii),flag_vec(ii)]=B_prima(B0,B_tilde,B_tilde,params);
    R_vec(ii)=1/beta*B_p_vec(ii)/B0;
end

figure
plot(B_vec,B_p_vec,'LineWidth',2); grid on; axis tight; hold on;
plot(B_vec,B_vec,'LineWidth',2,'LineStyle','--'); grid on; axis tight;
grid on; yvals=ylim; line([Bstar(B_tilde) Bstar(B_tilde)],[yvals(1) yvals(2)],'Color','k','LineWidth',2,'LineStyle','--');
grid on; yvals=ylim; line([Bast(B_tilde_p) Bast(B_tilde_p)],[yvals(1) yvals(2)],'Color','k','LineWidth',2,'LineStyle','--');
grid on; yvals=ylim; line([B_tilde B_tilde],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.5],'LineWidth',2,'LineStyle','-.');
grid on; yvals=ylim; line([Bhat(B_tilde,B_tilde) Bhat(B_tilde,B_tilde)],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.5],'LineWidth',2,'LineStyle','-.');
grid on; yvals=ylim; line([Bhat(B_tilde,B_tilde) Bhat(B_tilde,B_tilde)],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.5],'LineWidth',2,'LineStyle','-.');

figure
index1=(B_vec<Bstar(B_tilde));
index2=((B_vec>Bstar(B_tilde))&(B_vec<Bast(B_tilde)));
index3=(B_vec>Bast(B_tilde));
B_vec<Bast(B_tilde)
plot(B_vec(index1),R_vec(index1),'LineWidth',2,'Color',color_t0); axis tight; grid on; hold on;
plot(B_vec(index2),R_vec(index2),'LineWidth',2,'Color',color_t0); grid on;  
plot(B_vec(index3),R_vec(index3),'LineWidth',2,'Color',color_t0); grid on;  
title('$R(B)$','interpreter', 'latex'); 
% legend('$\mathcal{E}_{t}$','$\mathcal{E}_{t+1}$','Box','Off','interpreter', 'latex','Location','northwest','AutoUpdate','off')
grid on; yvals=ylim; line([Bstar(B_tilde) Bstar(B_tilde)],[yvals(1) yvals(2)],'Color','k','LineWidth',2,'LineStyle','--');
grid on; yvals=ylim; line([Bast(B_tilde) Bast(B_tilde)],[yvals(1) yvals(2)],'Color','k','LineWidth',2,'LineStyle','--');
grid on; yvals=ylim; line([B_tilde B_tilde],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.5],'LineWidth',2,'LineStyle','-.');
grid on; yvals=ylim; xvals=xlim; line([xvals(1) xvals(end)],[1/beta 1/beta],'Color',[0.3 0.3 0.3],'LineWidth',2,'LineStyle',':');
grid on; yvals=ylim; line([Bhat(B_tilde,B_tilde) Bhat(B_tilde,B_tilde)],[yvals(1) yvals(2)],'Color',[0.3 0.3 0.3],'LineWidth',2,'LineStyle','-.');
B_star_aux=Bstar(B_tilde);
set(gca,'XTick',B_star_aux,'XTickLabel',[]);
set(gca,'YTick',[yvals(1) yvals(2)],'YTickLabel',[]);
yTicks = get(gca,'ytick');
xTicks = get(gca, 'xtick');

ax = axis; %Get left most x-position
%text(ax(1)-0.01,yTicks(2),'$\mathcal{E}_t$', 'HorizontalAlignment','Right','interpreter', 'latex');
%text(ax(1)-0.15,yTicks(2),'$\mathcal{E}_{t+1}$', 'HorizontalAlignment','Right','interpreter', 'latex');
%text(ax(2),yTicks(1)-0.1,'$B_{t},B_{t+1}$', 'HorizontalAlignment','Right','interpreter', 'latex');
% text(ax(1)-0.01,1.45,'$y,i,\delta k$', 'HorizontalAlignment','Right','interpreter', 'latex');
text(ax(1)-0.01,yTicks(2),'$\beta^{-1}$', 'HorizontalAlignment','Right','interpreter', 'latex');
text(xTicks(1),yTicks(1)-0.01,'$B^{\star}$', 'HorizontalAlignment','Center','VerticalAlignment','Middle','interpreter', 'latex');
text(Bast(B_tilde),yTicks(1)-0.01,'$B^{\ast}(\tilde{B})$', 'HorizontalAlignment','Center','VerticalAlignment','Middle','interpreter', 'latex');
text(B_tilde,yTicks(1)-0.01,'$\tilde{B}$', 'HorizontalAlignment','Center','VerticalAlignment','Middle','interpreter', 'latex'); hold on;
if printit==1
    orient landscape;
   %  saveas(gcf,'F_SSRate','pdf');
   ax = gca;
    exportgraphics(ax,'F_SSRate.pdf','BackgroundColor','none');
end




%% Plot for R
% fix arbitrary debt level. 
B0=B_bar/5;
% B_vec=[linspace(B0-0.1,B0-0.01,10) linspace(B0,1+beta*B0+0.01,100) linspace(1+beta*B0+0.01,1+beta*B0+0.1,10)];
% Bp_vec=[linspace(B0-0.1,B0-0.01,10) linspace(B0,1+beta*B0+0.01,50) linspace(1+beta*B0+0.01,1+beta*B0+0.1,10)];
B_vec  = linspace(0.5,1+beta*B0+1,100)  ;
Bp_vec = linspace(0.5,1+beta*B0+1.1,50)   ;
B_p_mat = NaN(length(B_vec),length(Bp_vec));
R_mat   = NaN(length(B_vec),length(Bp_vec));
%Rg_mat=zeros(length(B_vec),length(Bp_vec))   ;
val_mat=zeros(length(B_vec),length(Bp_vec))  ;
flag_mat=zeros(length(B_vec),length(Bp_vec)) ;
%gamma_mat=zeros(length(B_vec),length(Bp_vec));
%rho_mat=zeros(length(B_vec),length(Bp_vec))  ;
%eta_mat=zeros(length(B_vec),length(Bp_vec))  ;
%one_mat=zeros(length(B_vec),length(Bp_vec))  ;


% Some Tests
%mu_1=mu_Z(B0,1+beta*B0);
%q_1 =q_Z(B0,1+beta*B0);
%mu_test = mu_Z(B0,B_vec); 
%ii=0; q_test=zeros(length(mu_test),1); q_test2=zeros(length(mu_test),1);
%for B_tilde=B_vec
%    ii=1+ii;
%    q_test(ii)  = q_Z(B0,B_tilde);
%    q_test2(ii) = (q_Z(B0,B_tilde)-1)*(B_tilde-B0);
%end

% Main Loop
ii=0;
%R_iter=1/beta;
for B_tilde=B_vec
    jj=0;
    ii=ii+1;
    for B_tilde_p=Bp_vec
        jj=jj+1;
        %if R_iter~=1/beta
        %    guess=R_iter*beta+0.01;
       % else
         %   guess=1+0.01;
        %end
        [B_p_mat(ii,jj),val_mat(ii,jj),flag_mat(ii,jj)]=B_prima(B0,B_tilde,B_tilde_p,params);
        R_mat(ii,jj)=1/beta*B_p_mat(ii,jj)/B0;
       % R_iter=R_mat(ii,jj);
    end
end

%% Key Figures - Phase Diagram
figure
surf(B_vec,Bp_vec,beta*R_mat');
title('$\beta \cdot R(B,\tilde{B},\tilde{B}^{''})$','Interpreter','latex');
xlabel('$\tilde{B}$','Interpreter','latex');
ylabel('$\tilde{B}^{''}$','Interpreter','latex'); axis tight;

figure
imagesc(B_vec,Bp_vec,R_mat');
title('$\beta \cdot R(B,\tilde{B},\tilde{B}^{''})$','Interpreter','latex');
xlabel('$\tilde{B}$','Interpreter','latex');
ylabel('$\tilde{B}^{''}$','Interpreter','latex'); axis tight;


% figure
% surf(B_vec,Bp_vec,eta_mat'); hold on; surf(Bp_vec,B_vec,gamma_mat); surf(Bp_vec,B_vec,eta_mat); surf(Bp_vec,B_vec,rho_mat);
% title('$\beta \cdot R(B,\tilde{B},\tilde{B}^{''})$','Interpreter','latex');
% xlabel('$\tilde{B}$','Interpreter','latex');
% ylabel('$\tilde{B}^{''}$','Interpreter','latex');

figure
surf(B_vec,Bp_vec,val_mat');
title('Numerical Error','Interpreter','latex');
xlabel('$\tilde{B}$','Interpreter','latex');
ylabel('$\tilde{B}^{''}$','Interpreter','latex');


figure
surf(B_vec,Bp_vec,flag_mat');
title('Flags','Interpreter','latex');
xlabel('$\tilde{B}$','Interpreter','latex');
ylabel('$\tilde{B}^{''}$','Interpreter','latex');

% surf(B_vec,Bp_vec,gamma_mat');
%     title('$\beta \cdot R(B,\tilde{B},\tilde{B}^{''})$','Interpreter','latex');
%     xlabel('$\tilde{B}$','Interpreter','latex');
%     ylabel('$\tilde{B}^{''}$','Interpreter','latex');
% 
% figure
%     surf(B_vec,Bp_vec,beta*Rg_mat');
%     title('$\beta \cdot R(B,\tilde{B},\tilde{B}^{''})$','Interpreter','latex','FontSize',fontsize);
%     xlabel('$\tilde{B}$','Interpreter','latex','FontSize',fontsize);
%     ylabel('$\tilde{B}^{''}$','Interpreter','latex','FontSize',fontsize); axis tight;
%     orient landscape;
%     saveas(gcf, 'fig_EqR', 'pdf')
% 
% 
% if plotit==1
%     figure
%     subplot(2,2,1)
%     surf(B_vec,Bp_vec,eta_mat');
%     title('$\beta \cdot R(B,\tilde{B},\tilde{B}^{''})$','Interpreter','latex');
%     xlabel('$\tilde{B}$','Interpreter','latex');
%     ylabel('$\tilde{B}^{''}$','Interpreter','latex');
%     
%     subplot(2,2,2)
%     surf(B_vec,Bp_vec,eta_mat');
%     title('$\beta \cdot R(B,\tilde{B},\tilde{B}^{''})$','Interpreter','latex');
%     xlabel('$\tilde{B}$','Interpreter','latex');
%     ylabel('$\tilde{B}^{''}$','Interpreter','latex');
%     
%     subplot(2,2,3)
%     surf(B_vec,Bp_vec,eta_mat');
%     title('$\beta \cdot R(B,\tilde{B},\tilde{B}^{''})$','Interpreter','latex');
%     xlabel('$\tilde{B}$','Interpreter','latex');
%     ylabel('$\tilde{B}^{''}$','Interpreter','latex');
%     
%     subplot(2,2,4)
%     surf(B_vec,Bp_vec,eta_mat');
%     title('$\beta \cdot R(B,\tilde{B},\tilde{B}^{''})$','Interpreter','latex');
%     xlabel('$\tilde{B}$','Interpreter','latex');
%     ylabel('$\tilde{B}^{''}$','Interpreter','latex');
% 
% end

%% Using B' Approach
% I'm using the following equation:
% (Bp - B * (1-beta*Bp)/(1-beta*B)*Pi)
% 

% params.B_p_res=B_p_res;
% 
% 
% Main plots
% range_test=[1.5 3]; B_tilde=2; B_tilde_p=2;
% figure
% fplot(@(x) q_E(x,B_tilde),range_test); title('Borrowing Marginal Price');
% 
% figure
% fplot(@(x) q_B(x,B_tilde),range_test); title('Borrowing Marginal Price');
% 
% figure
% fplot(@(x) Q_Z(x,B_tilde_p),range_test); title('Average Price');
% 
% figure
% fplot(@(x) Q_Z(x,B_tilde_p)./Q_Z(B0,B_tilde),range_test); title('Average Price Ratio');
% 
% figure
% fplot(@(x) Pi(2,x,B_tilde,B_tilde_p),range_test); title('Marginal Inflation');
% 
% figure
% fplot(@(x) -Euler_rhs(B0,x,B_tilde,B_tilde_p),range_test,'LineWidth',2); title(' Euler RHS (must be dec)');
% hold on;
% plot(range_test,range_test,'LineWidth',2);
% grid on; axis tight;
% 
% figure
% fplot(@(x) B_p_res(B0,x,B_tilde,B_tilde_p),range_test,'LineWidth',2); title('Euler RHS residual Euler Equation');
% hold on;
% plot(range_test,0*range_test,'LineWidth',2);
% grid on; axis tight;

%% Main Simulation
B_tilde_rho= 0.915 ; % mean reversal rate
B_tilde_init=0.0 ;
B_tilde_ss  = 0.1*B_bar        ;
B_star      = Bstar(B_tilde_ss);
B0          = B_star/2         ; 

% Example 1: Violent Crunch
T_pre      = 10  ; index_pre=1:T_pre;
T_crunch   = 20  ; index_crunch=T_pre+1:T_pre+T_crunch;
T_rec      = 60  ; index_rec=T_pre+T_crunch+1:T_pre+T_crunch+T_rec+1; 
T_post = 5; index_post = T_pre + T_crunch + T_rec + 1 : T_pre + T_crunch + T_rec + T_post;

% Initial Values
B_tilde_t   = [ones(1,T_pre)*B_tilde_ss zeros(1,T_crunch+T_rec)+B_tilde_init ones(1,T_post)*B_tilde_ss]';
B_tilde_v   = B_tilde_t;

% Example #1: Violent Transitions
% T = length(B_tilde_t);
sim_violent=sim_mod(B_tilde_t,params,B0);

% Example 2: Smooth Transition  
% B_tilde_t   = [B_tilde_ss*ones(1,T_pre) zeros(1,T_crunch) B_tilde_ss-B_tilde_rho.^(1:T_rec)*B_tilde_ss ones(1,T_post)*B_tilde_ss]';
B_tilde_t   = [ones(1,T_pre)*B_tilde_ss zeros(1,T_crunch)+B_tilde_init B_tilde_ss-B_tilde_rho.^(1:T_rec)*(B_tilde_ss-B_tilde_init) ones(1,T_post)*B_tilde_ss]';
B_tilde_s   = B_tilde_t;
T= length(B_tilde_t);
sim_smooth=sim_mod(B_tilde_t,params,B0);

% Plot Specs
linethick=3;
varlist={'B_tilde_t','B_t','R_t','C_s_t','S_w_t','mu_t','q_t','X_w_t','Y_t'};
vartaglist={'$\tilde{B}_t$','$B_t$','$R_t$','$C^s_t$','$S^w_t$','$\mu_t$','$q_t$','$X^w_t$','$Y_t$'};
printlist={'fig_tildeB','fig_B','fig_R','fig_C','fig_S','fig_mu','fig_q','fig_X','fig_Y'};
for vv=1:numel(varlist)
  %  eval(['var_v=sim_violent.' varlist{vv} ';']);
    eval(['var_s=sim_smooth.'  varlist{vv} ';']);
    vartag=vartaglist{vv};
    figure;
  %  plot(1:T,var_v,'LineWidth',linethick,'Color',color_t0); xlabel('Time','Interpreter','latex','FontSize',fontsize); 
    title(vartag,'Interpreter','latex','FontSize',fontsize); hold on; grid on; 
    % axis tight;
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

%% Planner Solutions
B_ss_violent=sim_violent.B_t(end);
theta_violent=1-(1-beta)*B_ss_violent;
sim_p_violent=sim_planner(B_tilde_v,theta_violent,params,B0);
color_t1=[0.6 0.4 0.4];

% Main Simulation
printlist={'fig2_tildeB','fig2_B','fig2_R','fig2_C','fig2_S','fig2_mu','fig2_q','fig2_X','fig2_Y'};
for vv=1:numel(varlist)
    eval(['var_v_ce=sim_violent.' varlist{vv} ';']);
    eval(['var_v_rp=sim_p_violent.'  varlist{vv} ';']);
    vartag=vartaglist{vv};
    figure;
    title(vartag,'Interpreter','latex','FontSize',fontsize); hold on; grid on; 
    plot(index_pre,var_v_ce(index_pre),'LineWidth',linethick,'LineStyle','-','Color',color_t0); hold on; 
    plot(index_pre,var_v_rp(index_pre),'LineWidth',linethick,'LineStyle','--','Color',color_t1); hold on; 
    scatter(index_pre(end),var_v_ce(index_pre(end)),'MarkerFaceColor','w','MarkerEdgeColor',color_t0); 
    plot(index_crunch-1,var_v_ce(index_crunch),'LineWidth',linethick,'LineStyle','-','Color',color_t0);
    plot(index_crunch-1,var_v_rp(index_crunch),'LineWidth',linethick,'LineStyle','--','Color',color_t1);
    scatter(index_crunch(end)-1,var_v_ce(index_crunch(end)),'MarkerFaceColor','w','MarkerEdgeColor',color_t0); 
    scatter(index_crunch(end)-1,var_v_ce(index_crunch(end)+1),'filled','MarkerFaceColor',color_t0,'MarkerEdgeColor',color_t0); 
    plot(index_rec-2,var_v_ce(index_rec),'LineWidth',linethick,'LineStyle','-','Color',color_t0);
    plot(index_rec-2,var_v_rp(index_rec),'LineWidth',linethick,'LineStyle','--','Color',color_t1);
    scatter(index_rec(end)-2,var_v_ce(index_rec(end)),'MarkerFaceColor','w','MarkerEdgeColor',color_t0); 
    scatter(index_rec(end)-2,var_v_ce(index_rec(end)),'MarkerFaceColor',color_t0,'MarkerEdgeColor',color_t0);
    plot(index_post-3,var_v_ce(index_post),'LineWidth',linethick,'LineStyle','-','Color',color_t0); grid on; % axis tight;
    plot(index_post-3,var_v_rp(index_post),'LineWidth',linethick,'LineStyle','--','Color',color_t1);
    scatter(index_post(end)-3,var_v_ce(index_post(end)),'MarkerFaceColor','w','MarkerEdgeColor',color_t0); 
    scatter(index_post(end)-3,var_v_ce(index_post(end)),'MarkerFaceColor',color_t0,'MarkerEdgeColor',color_t0);
    scatter(index_post(end)-3,var_v_ce(index_post(end)),'>','filled','MarkerFaceColor',color_t0);
    ybound=ylim;
    line([index_pre(end) index_pre(end)],ybound,'LineStyle',':','Color',[0.3 0.3 0.3]); 
    line([index_crunch(end)-1 index_crunch(end)-1],ybound,'LineStyle',':','Color',[0.3 0.3 0.3]); 
    line([index_rec(end)-2 index_rec(end)-2],ybound,'LineStyle',':','Color',[0.3 0.3 0.3]); 
    line([index_post(end)-3 index_post(end)-3],ybound,'LineStyle',':','Color',[0.3 0.3 0.3]);
    xlabel('Time','interpreter','latex');

    xlabel('Time','Interpreter','latex','FontSize',fontsize); title(vartag,'Interpreter','latex','FontSize',fontsize); hold on; 
    % axis tight;
   
    % grid on;
    if vv==1
        legend('Violent Transition (CE)','Violent Transition (PE)','Interpreter','Latex','FontSize',fontsize,'Location','NorthWest','Box','off');
    end
    orient landscape;
    eval(['saveas(gcf,''' printlist{vv} ''',''pdf'');']);
end

%% Smooth Transition
B_ss_smooth=sim_smooth.B_t(end);
theta_smooth=1-(1-beta)*B_ss_smooth;
sim_p_smooth=sim_planner(B_tilde_s,theta_smooth,params,B0);
color_t1=[0.6 0.4 0.4];
color_t1=[0.6 0.4 0.4];
% Main Simulation
printlist={'fig3_tildeB','fig3_B','fig3_R','fig3_C','fig3_S','fig3_mu','fig3_q','fig3_X','fig3_Y'};
for vv=1:numel(varlist)
    eval(['var_s_ce=sim_smooth.' varlist{vv} ';']);
    eval(['var_s_rp=sim_p_smooth.'  varlist{vv} ';']);
    vartag=vartaglist{vv};
    figure;
    plot(index_pre,var_s_ce(index_pre),'LineWidth',linethick,'LineStyle','-','Color',color_t0); hold on; 
    title(vartag,'Interpreter','latex','FontSize',fontsize); hold on; grid on; 
    plot(index_crunch-1,var_s_rp(index_crunch),'LineWidth',linethick,'LineStyle','--','Color',color_t1);
    % plot(index_pre,var_s_rp(index_pre),'LineWidth',linethick,'LineStyle','-','Color',color_t1); 
    scatter(index_pre(end),var_s_ce(index_pre(end)),'MarkerFaceColor','w','MarkerEdgeColor',color_t0); 
    plot(index_crunch-1,var_s_ce(index_crunch),'LineWidth',linethick,'LineStyle','-','Color',color_t0);
    plot(index_crunch-1,var_s_rp(index_crunch),'LineWidth',linethick,'LineStyle','--','Color',color_t1);
    scatter(index_crunch(end)-1,var_s_ce(index_crunch(end)),'MarkerFaceColor','w','MarkerEdgeColor',color_t0); 
    scatter(index_crunch(end)-1,var_s_ce(index_crunch(end)+1),'filled','MarkerFaceColor',color_t0,'MarkerEdgeColor',color_t0); 
    plot(index_rec-2,var_s_ce(index_rec),'LineWidth',linethick,'LineStyle','-','Color',color_t0);
    plot(index_rec-2,var_s_rp(index_rec),'LineWidth',linethick,'LineStyle','--','Color',color_t1);
    scatter(index_rec(end)-2,var_s_ce(index_rec(end)),'MarkerFaceColor','w','MarkerEdgeColor',color_t0); 
    scatter(index_rec(end)-2,var_s_ce(index_rec(end)),'MarkerFaceColor',color_t0,'MarkerEdgeColor',color_t0);
    plot(index_post-3,var_s_ce(index_post),'LineWidth',linethick,'LineStyle','-','Color',color_t0); grid on;%  axis tight;
    plot(index_post-3,var_s_rp(index_post),'LineWidth',linethick,'LineStyle','--','Color',color_t1);
    scatter(index_post(end)-3,var_s_ce(index_post(end)),'MarkerFaceColor','w','MarkerEdgeColor',color_t0); 
    scatter(index_post(end)-3,var_s_ce(index_post(end)),'MarkerFaceColor',color_t0,'MarkerEdgeColor',color_t0);
    scatter(index_post(end)-3,var_s_ce(index_post(end)),'>','filled','MarkerFaceColor',color_t0);
    ybound=ylim;
    line([index_pre(end) index_pre(end)],ybound,'LineStyle',':','Color',[0.3 0.3 0.3]); 
    line([index_crunch(end)-1 index_crunch(end)-1],ybound,'LineStyle',':','Color',[0.3 0.3 0.3]); 
    line([index_rec(end)-2 index_rec(end)-2],ybound,'LineStyle',':','Color',[0.3 0.3 0.3]); 
    line([index_post(end)-3 index_post(end)-3],ybound,'LineStyle',':','Color',[0.3 0.3 0.3]);
    xlabel('Time','interpreter','latex');
    % grid on;
    if vv==1
        legend('Competitive Equilbrium','Planner Solution','Interpreter','Latex','FontSize',fontsize,'Location','NorthWest','Box','off');
    end
    orient landscape;
    eval(['saveas(gcf,''' printlist{vv} ''',''pdf'');']);
end

%T=length(B_tilde_t);
%B_t=B0*ones(T,1);
%Rg_t=zeros(T,1); valg_t=zeros(T,1); flagg_t=zeros(T,1); 
%for tt=1:T-1
%    B0=B_t(tt); B_tilde=B_tilde_t(tt); B_tilde_p=B_tilde_t(tt+1);
%   [Rg_t(tt),valg_t(tt),flagg_t(tt)]=Rg(B0,B_tilde,B_tilde_p,params);
%   B_t(tt+1)=beta*Rg_t(tt)*B_t(tt); 
%end    

%T=length(B_tilde_t);
%B_t=B0*ones(T,1);
%Rg_t=zeros(T,1); valg_t=zeros(T,1); flagg_t=zeros(T,1); 

%for tt=1:T-1
%    B0=B_t(tt); B_tilde=B_tilde_t(tt); B_tilde_p=B_tilde_t(tt+1);
%    [Rg_t(tt),valg_t(tt),flagg_t(tt)]=Rg(B0,B_tilde,B_tilde_p,params);
%    B_t(tt+1)=b

% eta*Rg_t(tt)*B_t(tt); 
%end    


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

%% Key Functions
function A_out=CalA(mu,delta) 
    if mu==1
        A_out=10e-16;
    elseif mu==0
        A_out=delta;
    else
        A_out=delta/(1-delta)*(1-mu)./mu.*log((1-mu.*delta)./(1-mu));
    end
end

function q_out=q_mu(mu,delta) 
    T=length(mu); q_out=zeros(T,1);
    for tt=1:T
        if mu(tt)==1
            q_out(tt)=10e16;
        elseif mu(tt)==0
            q_out(tt)=1/delta;
        elseif mu(tt)<10e-10
            q_out(tt)=delta;
        else
            q_out(tt)=CalA(mu(tt),delta).^(-1);
        end
    end
end

function [B_p,fval,flag] = B_prima(B,B_tilde,B_tilde_p,params)
    Bstar= params.Bstar;
    Bast = params.Bast;
    B_p_res=params.B_p_res;
  %  if B_tilde<=B_tilde_p
        if B<=Bstar(B_tilde)
            B_p=B; fval=0; flag=1;
        elseif B>=Bast(B_tilde_p)
            B_p=B; fval=0; flag=1;
        else
            if B>Bstar(B_tilde)
                [B_p,fval,flag]=fsolve(@(x) B_p_res(B,x,B_tilde,B_tilde_p),B_tilde_p-0.1);
                B_p=max(B_p,Bstar(B_tilde));
            else
                %if B>Bast(B_tilde,B_tilde_p)
                [B_p,fval,flag]=fsolve(@(x) B_p_res(B,x,B_tilde,B_tilde_p),Bstar(B_tilde)-0.1);
                B_p=max(B_p,Bstar(B_tilde));
            %elseif 
             %   [B_p,fval,flag]=fsolve(@(x) B_p_res(B,x,B_tilde,B_tilde_p),B_tilde_p-0.2);
                 if Bstar(B_tilde_p)<B
                   B_p=max(B_p,Bstar(B_tilde_p));
                 end
            end
        end
%     elseif B_tilde>B_tilde_p
%         if B>=Bstar(B_tilde)
%             B_p=B; fval=0; flag=1;
%         elseif B<=Bast(B_tilde_p)
%             B_p=B_tilde_p; fval=0; flag=1;
%         else
%             [B_p,fval,flag]=fsolve(@(x) B_p_res(B,x,B_tilde,B_tilde_p),Bstar(B_tilde+0.1));
%             B_p=max(B_p,Bstar(B_tilde));
%         end
%      elseif B_tilde==B_tilde_p
%          B_p=B; fval=0; flag=0;
   % else
    %    B_p=NaN; fval=NaN; flag=NaN;
    % end
end

function sim=sim_mod(B_tilde_t,params,B0)
    % Begin main loop
    y_check=params.y_check;
    income_check=params.income_check;
    C_s=params.C_s;
    S_w=params.S_w;
    mu_Z=params.mu_Z;
    q_Z=params.q_Z;
    X_w=params.X_w;
    Y_Z=params.Y_Z;
    beta=params.beta; 
    delta=params.delta;

    % Begin main simulation
    T=length(B_tilde_t);
    B_t=B0*ones(T,1);
    Rg_t=ones(T,1); valg_t=zeros(T,1); flagg_t=zeros(T,1); 
    for tt=1:T-1
        B0=B_t(tt); B_tilde=B_tilde_t(tt); B_tilde_p=B_tilde_t(tt+1); % Update State Variable
        [B_t(tt+1),valg_t(tt),flagg_t(tt)]=B_prima(B0,B_tilde,B_tilde_p,params);
        Rg_t(tt)=B_t(tt+1)/B_t(tt); 
    end   
    C_s_t = C_s(B_t);
    S_w_t = S_w(B_t,B_tilde_t);
    mu_t  = mu_Z(B_t,B_tilde_t);
    q_t   = q_mu(mu_t,delta);
    X_w_t  = mu_t./q_t;
    Y_t   = Y_Z(B_t,B_tilde_t);
    y_check_t=y_check(B_t,B_tilde_t) ;
    income_check_t=income_check(B_t,B_tilde_t);

    % Probably we want to indicate the regions too

    % Saver Variables
    sim.B_t=B_t;
    sim.B_tilde_t=B_tilde_t;
    sim.R_t=1/beta*Rg_t;
    sim.C_s_t=C_s_t;
    sim.S_w_t=S_w_t;
    sim.mu_t=mu_t;
    sim.q_t=q_t;
    sim.X_w_t=X_w_t;
    sim.Y_t=Y_t;    
    sim.y_check_t=y_check_t;
    sim.income_check_t=income_check_t;
end

function sim=sim_planner(B_tilde_t,theta,params,B0)
    % Begin main loop
    y_check=params.y_check;
    income_check=params.income_check;
    C_s=params.C_s;
    S_w=params.S_w;
    mu_Z=params.mu_Z;
 %   q_Z=params.q_Z;
 %   X_w=params.X_w;
    Y_Z=params.Y_Z;
    P=params.P;
    B_Ram=params.B_Ram;

    % Parameters
    beta=params.beta; 
    delta=params.delta;

    % Begin main simulation
    T=length(B_tilde_t);
    B_t=B0*ones(T,1);
    Rg_t=ones(T,1); valg_t=zeros(T,1); flagg_t=zeros(T,1); 
    for tt=1:T
        if tt==1
            B0=B_t(tt); B_tilde=B_tilde_t(tt); B_tilde_p=B_tilde_t(tt+1); % Update State Variable
            [B_t(tt+1),valg_t(tt),flagg_t(tt)]=B_prima(B0,B_tilde,B_tilde_p,params);
            Rg_t(tt)=B_t(tt+1)/B_t(tt); 
        else
             B_tilde=B_tilde_t(tt); 
            B_t(tt)=B_Ram(B_tilde,theta);
        end
    end   
    C_s_t = C_s(B_t);
    S_w_t = S_w(B_t,B_tilde_t);
    mu_t  = mu_Z(B_t,B_tilde_t);
    q_t   = q_mu(mu_t,delta);
    X_w_t  = mu_t./q_t;
    Y_t   = Y_Z(B_t,B_tilde_t);
    y_check_t=y_check(B_t,B_tilde_t) ;
    income_check_t=income_check(B_t,B_tilde_t);

    % Probably we want to indicate the regions too

    % Saver Variables
    sim.B_t=B_t;
    sim.B_tilde_t=B_tilde_t;
    sim.R_t=1/beta*Rg_t;
    sim.C_s_t=C_s_t;
    sim.S_w_t=S_w_t;
    sim.mu_t=mu_t;
    sim.q_t=q_t;
    sim.X_w_t=X_w_t;
    sim.Y_t=Y_t;    
    sim.y_check_t=y_check_t;
    sim.income_check_t=income_check_t;
end

