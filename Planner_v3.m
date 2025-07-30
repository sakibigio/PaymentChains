%% Payment Crises: Code to Compute Equilibrium Solution 
% (c) Saki Bigio
clear; close all;

% Code Specs
plotit=0; % 1 if plot needed
printit=1;

%% Parameter Set
delta=0.9; % Value of discount in production
beta =0.95; % Discount factor
theta=0.75; % Pareto Weight
B_bar=1/(1-beta); 
fontsize=16     ;

%% Equilibrium Functions
% TFP Function
% CalA_ex=@(mu,delta) delta/(1-delta)*(1-mu)./mu.*log((1-mu.*delta)./(1-mu));
if plotit==1
    figure
    fplot(@(x) CalA_ex(x,0.9),[0.01 0.99],'LineWidth',3); xlabel('$\mu$','Interpreter','latex','FontSize',fontsize);
    title('$\mathcal{A}(\mu)$','interpreter','latex','FontSize',fontsize); grid on; axis tight; hold on;
    fplot(@(x) CalA_ex(x,0.5),[0.01 0.99],'LineWidth',3); xlabel('$\mu$','Interpreter','latex','FontSize',fontsize);
    legend('High $\delta$','Low $\delta$','Interpreter','Latex')
    orient landscape;
    saveas(gcf, 'TFP_example', 'pdf')
end

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

% Planner Objective
Bstar = @(B_tilde) max([1/beta*(B_tilde-1) 0]);
P    = @(B,B_tilde) (1-theta)*log((1-beta)*B)+theta*log(S_w(B,B_tilde)+X_w(B,B_tilde));
options = optimoptions('fminunc','Algorithm','quasi-newton');
% B_RamLeft = @(B_tilde) fmincon(@(x) -P(x,B_tilde),B_tilde-0.1,Bstar(B_tilde),B_tilde);
B_RamLeft = @(B_tilde) fminunc(@(x) -P(x,B_tilde),Bstar(B_tilde)*0.5+B_tilde*0.5);

B_RamRight= @(B_tilde) fmincon(@(x) -P(x,B_tilde),B_tilde+0.1,B_tilde,1/(1-beta));
% B_Ram= @(B_tilde) B_Ramright= @(B_tilde) fminunc(@(x) -P(x,B_tilde),B_tilde-0.0001,options);;

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

% Critical Points

Bast  = @(B_tilde_p) B_tilde_p/(q_Z(B_tilde_p,B_tilde_p).^(-1)-(q_Z(B_tilde_p,B_tilde_p).^(-1)-1)*(1-beta)*B_tilde_p);
Bhat  = @(B_tilde,B_tilde_p) fsolve(@(x) Euler_t(x,B_tilde)-1./(1/Bstar(B_tilde)-(1-beta)),Bstar(B_tilde)+0.2); 

% Right Elasticity 
Ell_Aux   = @(mu) (mu/(1-mu))*((1-delta)/(1-delta*mu)*log((1-delta*mu)/(1-mu))^(-1)-1);
Ell       = @(B_tilde) Ell_Aux(mu_Z(B,B_tilde));
EllA      = @(B,B_tilde) Ell_Aux(mu_Z(B,B_tilde))/q_Z(B,B_tilde);
RamLeftFOC= @(B,B_tilde) (1-(1-beta)*B)/((1-beta)*B)-theta/(1-theta)*Q_Z(B,B_tilde)*((1-beta*EllA(B,B_tilde))/(1-beta));

%% Checks: all satisfied
% make sure C+S+X adds to Y
y_check=@(B,B_tilde) Y_Z(B,B_tilde)-(C_s(B)+S_w(B,B_tilde)+X_w(B,B_tilde));

%% Optimal Unconstrained
B_opt  = (1-theta)/(1-beta)    ;
B_vec  = linspace(0.1,8,2000)  ;
P_out  = NaN(length(B_vec),1)  ;
P_bar  = P_out;

% Unconstrained Value
bb=0;
for B=B_vec
        bb=bb+1;
        P_bar(bb)= P(B,B_bar);
end
%B_opt_sx=B_Ramleft(B_bar);
check_best=P(B_opt,B_bar)-max(P_bar);

% Worse Value
P_underbar  = P_out;
% Solve 
bb=0;
for B=B_vec
        bb=bb+1;
        P_underbar(bb)= P(B,0);
end
B_opt_x=B_RamRight(0);
check_worse=P(B_opt_x,0)-max(P_underbar);

%% Plot Planner Objective
figure;
bt_vec=[theta-0.1:0.005:0.999]; 
B_optvec=zeros(1,length(bt_vec));
type_vec=B_optvec;
Bstar_vec=type_vec;
ii=0;
for bt_i=1:length(bt_vec);
    ii=ii+1;
    B_tilde= B_bar*(1-bt_vec(bt_i));
    Bstar_vec(bt_i)=Bstar(B_tilde);
    bb=0;
    for B=B_vec
        bb=bb+1;
        P_out(bb)= P(B,B_tilde);
    end
    if P(B_RamLeft(B_tilde),B_tilde)>P(B_RamRight(B_tilde),B_tilde);
        B_plan=B_RamLeft(B_tilde);
    else
        B_plan=B_RamRight(B_tilde);
    end
    B_optvec(bt_i)=B_plan;
    
    % Analytic Solution
    [B_p,P_val,type]= P_star(B_tilde,Bstar,B_opt,B_opt_x,RamLeftFOC,P,B_RamLeft,beta,delta,theta)
    B_optvec(bt_i)  = B_p;
    type_vec(bt_i)  = type;

    plot(B_vec,P_out,'LineWidth',1,'Color',[0.3+ii/length(bt_vec)/2 0.3+ii/length(bt_vec)/2 0.3+ii/length(bt_vec)/2]); hold on;
    if type==1
        scatter(B_p,P_val,'ob'); drawnow;
    elseif type==2
        scatter(B_p,P_val,'dr'); drawnow;
    elseif type==4
        scatter(B_p,P_val,'<m'); drawnow;
    elseif type==3
        scatter(B_p,P_val,'>g'); drawnow;
    end
    yvals=ylim; xvals=xlim;
    % set(gca,'XTick',Bstar(B_tilde),'XTickLabel',[]);
    % set(gca,'YTick',[yvals(1) yvals(2)],'YTickLabel',[]);
   % yTicks = get(gca,'ytick');
   % xTicks = get(gca, 'xtick');
   % grid on; yvals=ylim; line([Bstar(B_tilde) Bstar(B_tilde)],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.5],'LineWidth',2,'LineStyle',':');
    % grid on; yvals=ylim; line([Bast(B_tilde) Bast(B_tilde)],[yvals(1) yvals(2)],'Color','k','LineWidth',2,'LineStyle','--');
   % grid on; yvals=ylim; line([B_tilde B_tilde],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.5],'LineWidth',2,'LineStyle',':');
    % grid on; yvals=ylim; line([Bhat(B_tilde,B_tilde) Bhat(B_tilde,B_tilde)],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.5],'LineWidth',2,'LineStyle','-.');
    % grid on; yvals=ylim; line([Bhat(B_tilde,B_tilde) Bhat(B_tilde,B_tilde)],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.5],'LineWidth',2,'LineStyle','-.');
   % grid on; yvals=ylim; line([B_plan B_plan],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.6],'LineWidth',2,'LineStyle','-');
   % grid on; yvals=ylim; line([B_opt B_opt],[yvals(1) yvals(2)],'Color',[0.6 0.5 0.8],'LineWidth',2,'LineStyle',':');
    
   % text(Bstar(B_tilde),yTicks(1)-0.05,'$B^{\star}$', 'HorizontalAlignment','Center','VerticalAlignment','top','interpreter', 'latex','FontSize',14);
   % text(B_tilde,yTicks(1)-0.05,'$\tilde{B}$', 'HorizontalAlignment','Center','VerticalAlignment','top','interpreter', 'latex','FontSize',14); hold on;
   % text(B_opt,yTicks(1)-0.05,'$B_{ss}$', 'HorizontalAlignment','Center','VerticalAlignment','top','interpreter', 'latex','FontSize',14); hold on;
    % text(b_bar,y
    % text(b_bar,yTicks(1)-0.3,'$\bar{B}$', 'HorizontalAlignment','Center','VerticalAlignment','top','interpreter', 'latex','FontSize',14); hold on;

end
plot(B_vec,P_bar,'LineWidth',2,'Color',[0.7 0.7 0.8],'LineStyle',':'); %  title('Planner'); 
hold on;
plot(B_vec,P_underbar,'LineWidth',2,'Color',[0.8 0.7 0.7],'LineStyle',':'); % title('Planner'); 
hold off;

% Plotting Planner's Optimal B
figure
index_sub=(type_vec~=2);
Btilde_vec=B_bar*(1-bt_vec);
plot(Btilde_vec,Bstar_vec,'LineWidth',2.2,'Color',[0.7 0.8 0.8],'LineStyle',':'); hold on;
plot(Btilde_vec,Btilde_vec,'LineWidth',2.2,'Color',[0.8 0.5 0.5],'LineStyle',':'); axis tight;
plot(Btilde_vec(type_vec~=2),B_optvec(type_vec~=2),'LineWidth',2,'Color',[0.3 0.3 0.4],'LineStyle','-'); 
legend('$B^{\star}(\tilde{B})$','$\tilde{B}$','$B^{p}(\tilde{B})$','interpreter','latex','AutoUpdate','Off','Box','Off','Location','southeast')
plot(Btilde_vec(type_vec==2),B_optvec(type_vec==2),'LineWidth',2,'Color',[0.3 0.3 0.4],'LineStyle','-'); 
yvals=ylim; xvals=xlim;
yTicks = get(gca,'ytick');
xTicks = get(gca, 'xtick');
line([xvals(1) xvals(2)],[B_opt B_opt],'Color',[0.5 0.5 0.6],'LineWidth',1,'LineStyle',':');
line([xvals(1) xvals(2)],[B_opt_x B_opt_x],'Color',[0.5 0.5 0.6],'LineWidth',1,'LineStyle',':');
text(xvals(1)-0.3,B_opt,'$B_{ss}$', 'HorizontalAlignment','Center','VerticalAlignment','top','interpreter', 'latex','FontSize',14); hold on;
text(xvals(1)-0.3,B_opt_x,'$\underline{B}^{p}$', 'HorizontalAlignment','Center','VerticalAlignment','top','interpreter', 'latex','FontSize',14); hold on;
yTicks = get(gca,'ytick');
xTicks = get(gca, 'xtick');
set(gca,'XTick',B_tilde,'XTickLabel',[]);
set(gca,'YTick',[yvals(1) yvals(2)],'YTickLabel',[]);
xlabel('$\tilde{B}$','interpreter','latex','FontSize',14)
scatter(Btilde_vec(sum(index_sub)),B_optvec(sum(index_sub)),'MarkerEdgeColor',[0.3 0.3 0.4],'MarkerFaceColor',[0.3 0.3 0.4]);
scatter(Btilde_vec(sum(index_sub)+1),B_optvec(sum(index_sub)+1),'MarkerEdgeColor',[0.3 0.3 0.4],'MarkerFaceColor',[0.3 0.3 0.4]);
if printit==1
    orient landscape;
    ax = gca;
    exportgraphics(ax,'Planner_B.eps','BackgroundColor','none');
    % saveas(gcf, 'TFP_example', 'eps')
end

%% Single Plot - Planner Solution - Complement Example
% single plot
B_tilde= B_bar*(1-0.92);

bb=0;
for B=B_vec
    bb=bb+1;
    P_out(bb)= P(B,B_tilde);
end
if P(B_RamLeft(B_tilde),B_tilde)>P(B_RamRight(B_tilde),B_tilde);
    B_plan=B_RamLeft(B_tilde);
else
    B_plan=B_RamRight(B_tilde);
end
B_optvec(bt_i)=B_plan;

figure;
plot(B_vec,P_out,'LineWidth',3,'Color','k');  hold on; % 
grid on; xvals=xlim; line([xvals(1) xvals(2)],[max(P_underbar) max(P_underbar)],'Color',[0.5 0.5 0.8],'LineWidth',2,'LineStyle',':');
plot(B_vec,P_bar,'LineWidth',2,'Color',[0.5 0.5 0.8],'LineStyle',':'); hold on;
plot(B_vec,P_underbar,'LineWidth',2,'Color',[0.8 0.5 0.5],'LineStyle',':'); hold off;
legend('$\mathcal{P}(B,\tilde{B}_l)$','$\mathcal{P}(B,\bar{B})$','$\mathcal{P}(B,0)$','interpreter','latex','AutoUpdate','Off','Box','Off','Location','southeast')
yvals=ylim; xvals=xlim;
% set(gca,'XTick',Bstar(B_tilde),'XTickLabel',[]);
% set(gca,'YTick',[yvals(1) yvals(2)],'YTickLabel',[]);
yTicks = get(gca,'ytick');
xTicks = get(gca, 'xtick');
set(gca,'XTick',B_tilde,'XTickLabel',[]);
set(gca,'YTick',[yvals(1) yvals(2)],'YTickLabel',[]);
grid on; yvals=ylim; line([Bstar(B_tilde) Bstar(B_tilde)],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.5],'LineWidth',2,'LineStyle',':');
% grid on; yvals=ylim; line([Bast(B_tilde) Bast(B_tilde)],[yvals(1) yvals(2)],'Color','k','LineWidth',2,'LineStyle','--');
grid on; yvals=ylim; line([B_tilde B_tilde],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.5],'LineWidth',2,'LineStyle',':');
% grid on; yvals=ylim; line([Bhat(B_tilde,B_tilde) Bhat(B_tilde,B_tilde)],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.5],'LineWidth',2,'LineStyle','-.');
% grid on; yvals=ylim; line([Bhat(B_tilde,B_tilde) Bhat(B_tilde,B_tilde)],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.5],'LineWidth',2,'LineStyle','-.');
grid on; yvals=ylim; line([B_plan B_plan],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.6],'LineWidth',2,'LineStyle','-');
grid on; yvals=ylim; line([B_opt B_opt],[yvals(1) yvals(2)],'Color',[0.6 0.5 0.8],'LineWidth',2,'LineStyle',':');
% grid on; xvals=xlim; line([xvals(1) xvals(2)],[max(P_underbar) max(P_underbar)],'Color',[0.8 0.5 0.5],'LineWidth',2,'LineStyle',':');
% text(Bstar(B_tilde),yTicks(1)-0.2,'$B^{\star}$', 'HorizontalAlignment','Center','VerticalAlignment','top','interpreter', 'latex','FontSize',14);
text(Bstar(B_tilde),yTicks(1)-0.1,'$B^{\star}(\tilde{B})$', 'HorizontalAlignment','Center','VerticalAlignment','top','interpreter', 'latex','FontSize',14); hold on
text(B_tilde,yTicks(1)-0.1,'$\tilde{B}$', 'HorizontalAlignment','Center','VerticalAlignment','top','interpreter', 'latex','FontSize',14); hold on;
text(B_opt,yTicks(1)-0.1,'$B_{ss}$', 'HorizontalAlignment','Center','VerticalAlignment','top','interpreter', 'latex','FontSize',14); hold on;
text(B_plan,yTicks(1)-0.1,'$B^{p}$', 'HorizontalAlignment','Center','VerticalAlignment','top','interpreter', 'latex','FontSize',14); hold on;
grid on;
if printit==1
    orient landscape;
    ax = gca;
    exportgraphics(ax,'Planner_subs.eps','BackgroundColor','none');
    % saveas(gcf, 'TFP_example', 'eps')
end

%% Single Plot - Planner Solution - Complement Example
figure
plot(bt_vec,B_optvec);
% single plot
B_tilde= B_bar*(1-0.85);

bb=0;
for B=B_vec
    bb=bb+1;
    P_out(bb)= P(B,B_tilde);
end
if P(B_RamLeft(B_tilde),B_tilde)>P(B_RamRight(B_tilde),B_tilde);
    B_plan=B_RamLeft(B_tilde);
else
    B_plan=B_RamRight(B_tilde);
end
B_optvec(bt_i)=B_plan;

figure;
plot(B_vec,P_out,'LineWidth',3,'Color','k'); hold on; % title('Planner'); 
grid on; xvals=xlim; line([xvals(1) xvals(2)],[max(P_out) max(P_out)],'Color',[0.5 0.5 0.8],'LineWidth',2,'LineStyle',':');
plot(B_vec,P_bar,'LineWidth',2,'Color',[0.5 0.5 0.8],'LineStyle',':');  hold on; % title('Planner');
plot(B_vec,P_underbar,'LineWidth',2,'Color',[0.8 0.5 0.5],'LineStyle',':');  hold off; % title('Planner');
legend('$\mathcal{P}(B,\tilde{B}_l)$','$\mathcal{P}(B,\bar{B})$','$\mathcal{P}(B,0)$','interpreter','latex','AutoUpdate','Off','Box','Off','Location','southeast')
yvals=ylim; xvals=xlim;
% set(gca,'XTick',Bstar(B_tilde),'XTickLabel',[]);
% set(gca,'YTick',[yvals(1) yvals(2)],'YTickLabel',[]);
yTicks = get(gca,'ytick');
xTicks = get(gca, 'xtick');
set(gca,'XTick',B_tilde,'XTickLabel',[]);
set(gca,'YTick',[yvals(1) yvals(2)],'YTickLabel',[]);
grid on; yvals=ylim; line([Bstar(B_tilde) Bstar(B_tilde)],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.5],'LineWidth',2,'LineStyle',':');
% grid on; yvals=ylim; line([Bast(B_tilde) Bast(B_tilde)],[yvals(1) yvals(2)],'Color','k','LineWidth',2,'LineStyle','--');
grid on; yvals=ylim; line([B_tilde B_tilde],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.5],'LineWidth',2,'LineStyle',':');
% grid on; yvals=ylim; line([Bhat(B_tilde,B_tilde) Bhat(B_tilde,B_tilde)],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.5],'LineWidth',2,'LineStyle','-.');
% grid on; yvals=ylim; line([Bhat(B_tilde,B_tilde) Bhat(B_tilde,B_tilde)],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.5],'LineWidth',2,'LineStyle','-.');
grid on; yvals=ylim; line([B_plan B_plan],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.6],'LineWidth',2,'LineStyle','-');
grid on; yvals=ylim; line([B_opt B_opt],[yvals(1) yvals(2)],'Color',[0.6 0.5 0.8],'LineWidth',2,'LineStyle',':');
% text(Bstar(B_tilde),yTicks(1)-0.2,'$B^{\star}$', 'HorizontalAlignment','Center','VerticalAlignment','top','interpreter', 'latex','FontSize',14);
% text(Bstar(B_tilde),yTicks(1)-0.1,'$B^{\star}(\tilde{B})$', 'HorizontalAlignment','Center','VerticalAlignment','top','interpreter', 'latex','FontSize',14); hold on
text(B_tilde,yTicks(1)-0.1,'$\tilde{B}$', 'HorizontalAlignment','Center','VerticalAlignment','top','interpreter', 'latex','FontSize',14); hold on;
text(B_opt,yTicks(1)-0.1,'$B_{ss}$', 'HorizontalAlignment','Center','VerticalAlignment','top','interpreter', 'latex','FontSize',14); hold on;
text(B_plan,yTicks(1)-0.1,'$B^{p}$', 'HorizontalAlignment','Center','VerticalAlignment','top','interpreter', 'latex','FontSize',14); hold on;
grid on;
if printit==1
    orient landscape;
    ax = gca;
    exportgraphics(ax,'Planner_comp.eps','BackgroundColor','none');
    % saveas(gcf, 'TFP_example', 'eps')
end

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

%% Used Functions
function [B_p,P_val,type]=P_star(B_tilde,Bstar,B_opt,B_opt_x,RamLeftFOC,P,B_RamLeft,beta,delta,theta)
% type=1 uncontrained
% type=2 only chained
% type=3 constrained with some spot
% type=4 constrained with no chained
    if Bstar(B_tilde)>=B_opt
        B_p=B_opt;
        P_val=P(B_p,B_tilde);
        type=1;
    elseif RamLeftFOC(B_tilde,B_tilde)>0 
        B_p=B_opt_x;
        P_val=P(B_p,B_tilde);
        type=2;
    elseif  (1-(1-beta)*Bstar(B_tilde))/((1-beta)*Bstar(B_tilde))<theta/(1-theta)*(1-beta*delta)/(1-beta)
        B_p=Bstar(B_tilde);
        P_val=P(B_p,B_tilde);
        type=4;
        if P_val<=P(B_opt_x,B_tilde)
            B_p=B_opt_x;
            P_val=P(B_opt_x,B_tilde);
            type=2;
        end
    else
        B_p=B_RamLeft(B_tilde);
        P_val=P(B_p,B_tilde);
        type=3;
        if P_val<=P(B_opt_x,B_tilde)
            B_p=B_opt_x;
            P_val=P(B_opt_x,B_tilde);
            type=2;
        end
    end
end

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