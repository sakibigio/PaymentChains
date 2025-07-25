%% Payment Chains
% (c) 
% version with only consumption, no production. 
clear; close all;
tol=10e-6;
printinfo=0;
printit=0;

%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%
% Value functions - Partial Equilibrium
beta= 0.8   ; % 0.8
R   = 1/beta;
q   = 1.75  ; % 1.75
b_bar=1/(1-beta); % one unit of labor input, or annuity value: natural borrowing limit
B_tilde=0.4*b_bar; % 0.4
Bstar=1/beta*(B_tilde-1); % Threshold
Bast  = B_tilde/(q^(-1)-(q^(-1)-1)*(1-beta)*B_tilde);

% Suppor of values of debt
range=[-1 b_bar-10e-10]    ;
B_vec=linspace(range(1),range(end),1000); % grid space

%% Special Solutions
% Solution without payment networks: exogenous prices
val_high=@(B) log((1-(1-beta)*B))/(1-beta)  ; % no chained consumption
val_low=@(B) log((1-(1-beta)*B)/q)/(1-beta) ; % all consumption is at chained price q
C_low=@(B) (1-(1-beta)*B)/q; % annuity value of consumption

% figure
% fplot(@(B) val_high(B),range,'LineWidth',3,'LineStyle','-'); hold on; axis tight;
% fplot(@(B) val_low(B),range,'LineWidth',3,'LineStyle','-'); grid on;% fplot(@(B) val_high(B)+log(1/q)/(1-beta),range,'LineWidth',3,'LineStyle',':'); hold on; axis tight;



%% Value Function Iteration
% define functions - B_p is future debt, determining expenditures today
E     = @(B_p,B) 1-B+B_p/R                    ;
% knowing expenditures and credit line, you can determine optimal
% expenditures in spot and chained, and determine total consumption
S_w   = @(E,B,B_tilde) min(max(B_tilde-B,0),E);
X_w   = @(E,B,B_tilde) (E-S_w(E,B,B_tilde))/q ;
C_w   = @(E,B,B_tilde) X_w(E,B,B_tilde)+S_w(E,B,B_tilde)          ;

% Main Iteration 
% initialize values
V0=val_low(B_vec); V_out=V0; C_out=0*V_out; S_out=0*V_out; X_out=0*V_out;; Bp_out=0*V_out;
% reset tolerance to enter loop
val_gap=2*tol;
while val_gap>tol
    if printinfo
        plot(B_vec,V0,'LineWidth',1,'LineStyle',':'); drawnow;
        disp(['Val at iter:' num2str(val_gap)]);
    end
    bb=0;
    for B=B_vec
        bb=bb+1;
        % maximize value function based on cubic spline
        Tv=@(B_p,B) log(C_w(E(B_p,B),B,B_tilde))+beta*interp1(B_vec,V0,B_p,'cubic');
        % compute all possible values
        Tv_iter=Tv(B_vec,B);
        index=(imag(Tv_iter)==0);
        % evaluate at optimum to obtain solution
        [V_out(bb),I]=max(Tv_iter(index));
        Bp_aux=B_vec(index);
        % update solutions
        C_out(bb)=C_w(E(Bp_aux(I),B),B,B_tilde);
        S_out(bb)=S_w(E(Bp_aux(I),B),B,B_tilde);
        X_out(bb)=X_w(E(Bp_aux(I),B),B,B_tilde);
        Bp_out(bb)=Bp_aux(I);
    end
    val_gap=abs(max(V_out-V0));
    V0=V_out;
end
% Numerical Exit point
B_exit_num=min(B_vec(abs(V0-val_low(B_vec))<10e-10));

%% Computing Conjectured Value after which you do not exit
Tv_iter=Tv(B_vec,Bast);
[V_Bast,I]=max(Tv_iter(Tv_iter==real(Tv_iter)));
%B_exit=fsolve(@(b) log((1-b+Bast/R)/q)+beta*V_Bast-val_low(b),Bast);
% B_exit=fsolve(@(b) log(C_w(E(Bast,b),b,B_tilde))+beta*interp1(B_vec,V0,Bast,'cubic'),Bast);
% B_eUB=fsolve(@(b)  log((Bast-1/beta*b+1)/q)+beta*val_high(Bast)-val_low(b),Bast);

% Third Posibility - Value function that jumps to above threshold in one
% shot...
% B_exit=fsolve(@(b) log((1-(1/beta*b-(Bast-0.1)))/q)+beta*interp1(B_vec,V0,Bast-0.1)-val_low(b),Bast)

% Plot Solutions - Worker Problem
figure
% Colors
color_t1=[0.7 0.7 0.7];
color_t0=[0.0 0.0 0.6];
fplot(@(B) val_high(B),range,'LineWidth',1.5,'LineStyle','-','Color',color_t0); hold on; axis tight;
fplot(@(B) val_low(B),range,'LineWidth',1.5,'LineStyle','-','Color',color_t1); grid on;
% fplot(@(B) val_high(B)+log(1/q)/(1-beta),range,'LineWidth',3,'LineStyle',':','Color',color_t1); hold on; axis tight;
plot(B_vec,V0,'LineWidth',2,'LineStyle','-','Color','k'); drawnow;
legend('$\bar{V}$','$\underbar{V}$','$V_ss(B,\tilde{B})$','Interpreter','Latex','AutoUpdate','off','Box','Off')
yvals=ylim; xvals=xlim;
set(gca,'XTick',Bstar,'XTickLabel',[]);
set(gca,'YTick',[yvals(1) yvals(2)],'YTickLabel',[]);
yTicks = get(gca,'ytick');
xTicks = get(gca, 'xtick');
grid on; yvals=ylim; line([Bstar Bstar],[yvals(1) yvals(2)],'Color','k','LineWidth',1,'LineStyle','--');
grid on; yvals=ylim; line([Bast Bast],[yvals(1) yvals(2)],'Color','k','LineWidth',1,'LineStyle','--');
grid on; yvals=ylim; line([B_tilde B_tilde],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.5],'LineWidth',1,'LineStyle','-.'); grid minor;
% grid on; yvals=ylim; line([B_exit B_exit],[yvals(1) yvals(2)],'Color',[0.5 0.0 0.0],'LineWidth',1,'LineStyle','-.'); grid minor;
grid on; yvals=ylim; line([B_exit_num B_exit_num],[yvals(1) yvals(2)],'Color',[0.5 0.0 0.0],'LineWidth',1,'LineStyle','-.'); grid minor;

% grid on; yvals=ylim; line([Bhat(B_tilde,B_tilde) Bhat(B_tilde,B_tilde)],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.5],'LineWidth',1,'LineStyle',':');
text(xTicks(1),yTicks(1)-1,'$B^{\star}$', 'HorizontalAlignment','Center','VerticalAlignment','top','interpreter', 'latex','FontSize',14);
text(Bast,yTicks(1)-1,'$B^{\ast}(\tilde{B})$', 'HorizontalAlignment','Center','VerticalAlignment','top','interpreter', 'latex','FontSize',14);
text(B_tilde,yTicks(1)-1,'$\tilde{B}$', 'HorizontalAlignment','Center','VerticalAlignment','top','interpreter', 'latex','FontSize',14); hold on;
text(b_bar,yTicks(1)-1,'$\bar{B}$', 'HorizontalAlignment','Center','VerticalAlignment','top','interpreter', 'latex','FontSize',14); hold on;
text(B_exit_num,yTicks(1)-1,'$B^{h}$', 'HorizontalAlignment','Center','VerticalAlignment','top','interpreter', 'latex','FontSize',14); hold on;
printit=0; grid on;
if printit==1
    orient landscape;
    % saveas(gcf,'F_valuefunction','pdf');
    ax= gca;
    exportgraphics(ax,'F_valuefunction.pdf','BackgroundColor','none');
end

%% Consumption Policy Actual
figure
% Colors
color_t1=[0.7 0.7 0.7];
color_t0=[0.0 0.0 0.6];
%fplot(@(B) val_high(B),range,'LineWidth',1.5,'LineStyle','-','Color',color_t0); hold on; axis tight;
% fplot(@(B) val_low(B),range,'LineWidth',1.5,'LineStyle','-','Color',color_t1); grid on;
% fplot(@(B) val_high(B)+log(1/q)/(1-beta),range,'LineWidth',3,'LineStyle',':','Color',color_t1); hold on; axis tight;
% plot(B_vec,C_high(B_vec),'LineWidth',2,'LineStyle',':','Color','k'); 
plot(B_vec,C_out,'LineWidth',2,'LineStyle','-','Color','k'); drawnow; grid on; hold on;
plot(B_vec,C_low(B_vec),'LineWidth',2,'LineStyle',':','Color','k'); 
legend('$C(B;\tilde{B}=0)$','$C(B;\tilde{B}>0)$','Interpreter','Latex','AutoUpdate','off','Box','Off')
yvals=ylim; xvals=xlim;
set(gca,'XTick',Bstar,'XTickLabel',[]);
set(gca,'YTick',[yvals(1) yvals(2)],'YTickLabel',[]);
yTicks = get(gca,'ytick');
xTicks = get(gca, 'xtick');
grid on; yvals=ylim; line([Bstar Bstar],[yvals(1) yvals(2)],'Color','k','LineWidth',1,'LineStyle','--');
grid on; yvals=ylim; line([Bast Bast],[yvals(1) yvals(2)],'Color','k','LineWidth',1,'LineStyle','--');
grid on; yvals=ylim; line([B_tilde B_tilde],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.5],'LineWidth',1,'LineStyle','-.'); grid minor;
% grid on; yvals=ylim; line([B_exit B_exit],[yvals(1) yvals(2)],'Color',[0.5 0.0 0.0],'LineWidth',1,'LineStyle','-.'); grid minor;
grid on; yvals=ylim; line([B_exit_num B_exit_num],[yvals(1) yvals(2)],'Color',[0.5 0.0 0.0],'LineWidth',1,'LineStyle','-.'); grid minor;

% grid on; yvals=ylim; line([Bhat(B_tilde,B_tilde) Bhat(B_tilde,B_tilde)],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.5],'LineWidth',1,'LineStyle',':');
text(xTicks(1),yTicks(1)-0.1,'$B^{\star}$', 'HorizontalAlignment','Center','VerticalAlignment','bottom','interpreter', 'latex');
text(Bast,yTicks(1)-0.1,'$B^{\ast}(\tilde{B})$', 'HorizontalAlignment','Center','VerticalAlignment','bottom','interpreter', 'latex');
text(B_tilde,yTicks(1)-0.1,'$\tilde{B}$', 'HorizontalAlignment','Center','VerticalAlignment','bottom','interpreter', 'latex'); hold on;


%% Consumption Policy - Split
figure
% Colors
color_t1=[0.7 0.7 0.7];
color_t0=[0.0 0.0 0.6];
%fplot(@(B) val_high(B),range,'LineWidth',1.5,'LineStyle','-','Color',color_t0); hold on; axis tight;
% fplot(@(B) val_low(B),range,'LineWidth',1.5,'LineStyle','-','Color',color_t1); grid on;
% fplot(@(B) val_high(B)+log(1/q)/(1-beta),range,'LineWidth',3,'LineStyle',':','Color',color_t1); hold on; axis tight;
% plot(B_vec,C_high(B_vec),'LineWidth',2,'LineStyle',':','Color','k'); 
plot(B_vec,S_out,'LineWidth',2,'LineStyle','-','Color','k'); drawnow; grid on; hold on;
plot(B_vec,X_out,'LineWidth',2,'LineStyle','-','Color',[0.7 0.7 0.7]); drawnow; grid on; hold on;
plot(B_vec,C_low(B_vec),'LineWidth',2,'LineStyle',':','Color','k'); 
legend('$S^{w}(B)$','$X^{w}(B)$','$\underbar{C}(B)$','Interpreter','Latex','AutoUpdate','off','Box','Off')
yvals=ylim; xvals=xlim;
set(gca,'XTick',Bstar,'XTickLabel',[]);
set(gca,'YTick',[yvals(1) yvals(2)],'YTickLabel',[]);
yTicks = get(gca,'ytick');
xTicks = get(gca, 'xtick');
yvals=ylim; line([Bstar Bstar],[yvals(1) yvals(2)],'Color','k','LineWidth',1,'LineStyle','--');
yvals=ylim; line([Bast Bast],[yvals(1) yvals(2)],'Color','k','LineWidth',1,'LineStyle','--');
yvals=ylim; line([B_tilde B_tilde],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.5],'LineWidth',1,'LineStyle','-.'); grid minor;
% grid on; yvals=ylim; line([B_exit B_exit],[yvals(1) yvals(2)],'Color',[0.5 0.0 0.0],'LineWidth',1,'LineStyle','-.'); grid minor;
 yvals=ylim; line([B_exit_num B_exit_num],[yvals(1) yvals(2)],'Color',[0.5 0.0 0.0],'LineWidth',1,'LineStyle','-.'); grid minor;

% grid on; yvals=ylim; line([Bhat(B_tilde,B_tilde) Bhat(B_tilde,B_tilde)],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.5],'LineWidth',1,'LineStyle',':');
text(xTicks(1),yTicks(1)-0.205,'$B^{\star}$', 'HorizontalAlignment','Center','VerticalAlignment','top','interpreter', 'latex','FontSize',14);
text(Bast,yTicks(1)-0.025,'$B^{\ast}(\tilde{B})$', 'HorizontalAlignment','Center','VerticalAlignment','top','interpreter', 'latex','FontSize',14);
text(B_tilde,yTicks(1)-0.025,'$\tilde{B}$', 'HorizontalAlignment','Center','VerticalAlignment','top','interpreter', 'latex','FontSize',14); hold on;
text(b_bar,yTicks(1)-0.025,'$\bar{B}$', 'HorizontalAlignment','Center','VerticalAlignment','top','interpreter', 'latex','FontSize',14); hold on;
text(B_exit_num,yTicks(1)-0.025,'$B^{h}$', 'HorizontalAlignment','Center','VerticalAlignment','top','interpreter', 'latex','FontSize',14); hold on;
grid on;
if printit==1
    orient landscape;
    ax = gca;
    exportgraphics(ax,'F_cons.pdf','BackgroundColor','none');
%     saveas(gcf,'F_cons','pdf');
end


%% Debt Policy
figure
% Colors
color_t1=[0.7 0.7 0.7];
color_t0=[0.0 0.0 0.6];
%fplot(@(B) val_high(B),range,'LineWidth',1.5,'LineStyle','-','Color',color_t0); hold on; axis tight;
% fplot(@(B) val_low(B),range,'LineWidth',1.5,'LineStyle','-','Color',color_t1); grid on;
% fplot(@(B) val_high(B)+log(1/q)/(1-beta),range,'LineWidth',3,'LineStyle',':','Color',color_t1); hold on; axis tight;
plot(B_vec,Bp_out,'LineWidth',2,'LineStyle','-','Color','k'); grid on; drawnow; hold on;
plot(B_vec,B_vec,'LineWidth',2,'LineStyle',':','Color','k');
legend('$B''$','$B$','Interpreter','Latex','AutoUpdate','off','Box','Off','Location','northwest')
yvals=ylim; xvals=xlim;
set(gca,'XTick',Bstar,'XTickLabel',[]);
set(gca,'YTick',[yvals(1) yvals(2)],'YTickLabel',[]);
yTicks = get(gca,'ytick');
xTicks = get(gca, 'xtick');
line([xvals(1) xvals(2)],[Bstar Bstar],'Color','k','LineWidth',1,'LineStyle',':');
line([xvals(1) xvals(2)],[B_tilde B_tilde],'Color','k','LineWidth',1,'LineStyle',':');

yvals=ylim; line([Bstar Bstar],[yvals(1) yvals(2)],'Color','k','LineWidth',1,'LineStyle','--');
yvals=ylim; line([Bast Bast],[yvals(1) yvals(2)],'Color','k','LineWidth',1,'LineStyle','--');
yvals=ylim; line([B_tilde B_tilde],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.5],'LineWidth',1,'LineStyle','-.'); % grid minor;
yvals=ylim; line([B_exit_num B_exit_num],[yvals(1) yvals(2)],'Color',[0.5 0.0 0.0],'LineWidth',1,'LineStyle','-.'); % grid minor;
% grid on; yvals=ylim; line([Bhat(B_tilde,B_tilde) Bhat(B_tilde,B_tilde)],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.5],'LineWidth',1,'LineStyle',':');
text(xTicks(1),yTicks(1)-0.3,'$B^{\star}$', 'HorizontalAlignment','Center','VerticalAlignment','top','interpreter', 'latex','FontSize',14);
text(Bast,yTicks(1)-0.3,'$B^{\ast}(\tilde{B})$', 'HorizontalAlignment','Center','VerticalAlignment','top','interpreter', 'latex','FontSize',14);
text(B_tilde,yTicks(1)-0.3,'$\tilde{B}$', 'HorizontalAlignment','Center','VerticalAlignment','top','interpreter', 'latex','FontSize',14); hold on;
text(b_bar,yTicks(1)-0.3,'$\bar{B}$', 'HorizontalAlignment','Center','VerticalAlignment','top','interpreter', 'latex','FontSize',14); hold on;
text(B_exit_num,yTicks(1)-0.3,'$B^{h}$', 'HorizontalAlignment','Center','VerticalAlignment','top','interpreter', 'latex','FontSize',14); hold on;
if printit==1    
    orient landscape;
  %  set(gca,'LooseInset',get(gca,'TightInset'));
  %  saveas(gcf,'F_bprime','pdf');
    ax = gca;
    exportgraphics(ax,'F_bprime.pdf','BackgroundColor','none');
end

%% Backwards Solution Using Euler Equation
% Better done with a while loop
c_star=1+Bstar*(1/R-1)
c_path=c_star;
E_path=c_path;
Q_path=0*c_path+1;
B_path=Bstar;
B_aux=B_tilde-0.001;
while B_path(1)<B_tilde
    c_aux=(beta*R*q)^(-1)*c_path(1);
    B_aux=1/q+B_path(1)/(q*R)-c_aux+(q-1)*B_tilde/q;
    Q_aux=-(q-1)*(B_tilde-B_aux)/c_aux+q;
    B_check=(1+B_path(1)/R-Q_aux*c_aux)-B_aux;
    if B_aux<B_tilde
        c_path=[c_aux; c_path];    
        B_path=[B_aux; B_path];
        Q_path=[Q_aux; c_path];
    else
        uB_aux=1+B_path(1)/R-q*c_aux;
        c_path=[c_aux; c_path];    
        B_path=[B_aux; B_path];
        Q_path=[Q_aux; c_path];
        break
    end
end
iter=0;
c_aux=(beta*R*q)^(-1)*c_path(2);
while B_path(1)<B_exit_num
    B_aux=1+B_path(1)/R-q*c_aux;
    Q_aux=q;
    c_path=[c_aux; c_path];    
    B_path=[B_aux; B_path];
    Q_path=[Q_aux; c_path];
    if length(c_path)>20
        break
    end
    c_aux=(beta*R)^(-1)*c_path(1);
end
T_trans=length(c_path)-1;

figure
% Colors
color_t1=[0.7 0.7 0.7];
color_t0=[0.0 0.0 0.6];
%fplot(@(B) val_high(B),range,'LineWidth',1.5,'LineStyle','-','Color',color_t0); hold on; axis tight;
% fplot(@(B) val_low(B),range,'LineWidth',1.5,'LineStyle','-','Color',color_t1); grid on;
% fplot(@(B) val_high(B)+log(1/q)/(1-beta),range,'LineWidth',3,'LineStyle',':','Color',color_t1); hold on; axis tight;
% plot(B_vec,C_high(B_vec),'LineWidth',2,'LineStyle',':','Color','k'); 
plot(B_vec,C_out,'LineWidth',2,'LineStyle','-','Color','k'); drawnow; grid on; hold on;
legend('$C(B;\tilde{B}>0)$','Interpreter','Latex','AutoUpdate','off','Box','Off')
yvals=ylim; xvals=xlim;
set(gca,'XTick',Bstar,'XTickLabel',[]);
set(gca,'YTick',[yvals(1) yvals(2)],'YTickLabel',[]);
yTicks = get(gca,'ytick');
xTicks = get(gca, 'xtick');
grid on; yvals=ylim; line([Bstar Bstar],[yvals(1) yvals(2)],'Color','k','LineWidth',1,'LineStyle','--');
grid on; yvals=ylim; line([Bast Bast],[yvals(1) yvals(2)],'Color','k','LineWidth',1,'LineStyle','--');
grid on; yvals=ylim; line([B_tilde B_tilde],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.5],'LineWidth',1,'LineStyle','-.'); grid minor;
% grid on; yvals=ylim; line([B_exit B_exit],[yvals(1) yvals(2)],'Color',[0.5 0.0 0.0],'LineWidth',1,'LineStyle','-.'); grid minor;
grid on; yvals=ylim; line([B_exit_num B_exit_num],[yvals(1) yvals(2)],'Color',[0.5 0.0 0.0],'LineWidth',1,'LineStyle','-.'); grid minor;
for TT=(T_trans:-1:2)
    scatter(B_path(TT),c_path(TT),'r*');
end
% grid on; yvals=ylim; line([Bhat(B_tilde,B_tilde) Bhat(B_tilde,B_tilde)],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.5],'LineWidth',1,'LineStyle',':');
text(xTicks(1),yTicks(1)-0.1,'$B^{\star}$', 'HorizontalAlignment','Center','VerticalAlignment','bottom','interpreter', 'latex');
text(Bast,yTicks(1)-0.1,'$B^{\ast}(\tilde{B})$', 'HorizontalAlignment','Center','VerticalAlignment','bottom','interpreter', 'latex');
text(B_tilde,yTicks(1)-0.1,'$\tilde{B}$', 'HorizontalAlignment','Center','VerticalAlignment','bottom','interpreter', 'latex'); hold on;


%% Computing Backward Sets
% Compute Maximum T
cond=0; T=0; gamma=beta/q;
while cond==0
    T=T+1;
    left=0;
    for tt=0:T-1
        left=left+gamma^(tt)*(1-q^(-T+tt));
    end
    right=0;
    for tt=0:T
        right=right+gamma^(tt)*(1-q^(-T+tt-1));
    end
    if left<1&&right>1
        cond=1;
    end
%     if T>10;
%         break;
%     end
end

% alternative sum
cond=0; T=0;
scale_l=@(t) (1-gamma^(t))/(1-gamma)-q^(-t)*(1-beta^t)/(1-beta);
while cond==0
    T=T+1;
    left=(1-gamma^(T))/(1-gamma)-q^(-T)*(1-beta^T)/(1-beta);
    right=(1-gamma^(T+1))/(1-gamma)-q^(-(T+1))*(1-beta^(T+1))/(1-beta);
    if left<1&&right>=1
        cond=1;
    end
%     if T>1000;
%         break;
%     end
end

% Compute Backward Sets - entering mixed region
B_back_set=zeros(T+1,1);
B_l=@(T) Bstar+(B_tilde-Bstar)*scale_l(T);
for tt=1:T+1
    B_back_set(tt)=B_l(tt);
end


% Compute Backward Sets - chained region
N_tau=100;
B_back_l_set=zeros(N_tau,1);
B_back_r_set=zeros(N_tau,1);
consf_l=1-q^(-(T))*(B_tilde-Bstar);
consf_r=1-q^(-(T))*(B_tilde-Bstar);
%consf_l_test=1-q^(-T)*(1-(1-beta)*Bstar);
B_l2=@(t) beta^(t+1)*B_l(T)+(1-beta^(t+1))/(1-beta)*consf_l;
B_r2=@(t) beta^(t+1)*B_l(T+1)+(1-beta^(t+1))/(1-beta)*consf_r;
for tt=1:N_tau
    B_back_l_set(tt)=B_l2(tt);
    B_back_r_set(tt)=B_r2(tt);
end
% B_back_lim=1/(1-beta)*consf_l;
B_back_lim=1/(1-beta)*consf_r;

aux=0;
for tt=1:T
    aux=aux+beta^(tt-1)*log(q^(tt)*q^(-(T+1))*(B_tilde-Bstar));
end

V_d_bound=zeros(N_tau,1);
for tt=1:N_tau    
    V_d_bound(tt)=beta^(tt+1)*aux+(1-beta^(tt+1))*log(q^(-(T+1))*(B_tilde-Bstar))/(1-beta)+beta^(T+tt+1)*log((B_tilde-Bstar))/(1-beta);
end



figure('Name','Backward Sets','NumberTitle','off')
% Colors
color_t1=[0.7 0.7 0.7];
color_t0=[0.0 0.0 0.6];
fplot(@(B) val_high(B),range,'LineWidth',1.5,'LineStyle','-','Color',color_t0); hold on; axis tight;
fplot(@(B) val_low(B),range,'LineWidth',1.5,'LineStyle','-','Color',color_t1); grid on;
% fplot(@(B) val_high(B)+log(1/q)/(1-beta),range,'LineWidth',3,'LineStyle',':','Color',color_t1); hold on; axis tight;
plot(B_vec,V0,'LineWidth',2,'LineStyle','-','Color','k'); drawnow;
legend('$\bar{V}$','$\underbar{V}$','$V_ss(B,\tilde{B})$','Interpreter','Latex','AutoUpdate','off','Box','Off')
yvals=ylim; xvals=xlim;
set(gca,'XTick',Bstar,'XTickLabel',[]);
set(gca,'YTick',[yvals(1) yvals(2)],'YTickLabel',[]);
yTicks = get(gca,'ytick');
xTicks = get(gca, 'xtick');
grid on; yvals=ylim; line([Bstar Bstar],[yvals(1) yvals(2)],'Color','k','LineWidth',1,'LineStyle','--');
% Hi grid on; yvals=ylim; line([Bast Bast],[yvals(1) yvals(2)],'Color','k','LineWidth',1,'LineStyle','--');
grid on; yvals=ylim; line([B_tilde B_tilde],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.5],'LineWidth',1,'LineStyle','-.'); grid minor;
% grid on; yvals=ylim; line([B_exit B_exit],[yvals(1) yvals(2)],'Color',[0.5 0.0 0.0],'LineWidth',1,'LineStyle','-.'); grid minor;
grid on; yvals=ylim; line([B_exit_num B_exit_num],[yvals(1) yvals(2)],'Color',[0.5 0.0 0.0],'LineWidth',1,'LineStyle','-.'); grid minor;
%for tt=1:N_tau
text(xTicks(1),yTicks(1)-0.3,'$B^{\star}$', 'HorizontalAlignment','Center','VerticalAlignment','top','interpreter', 'latex','FontSize',14);
% text(Bast,yTicks(1)-0.3,'$B^{\ast}(\tilde{B})$', 'HorizontalAlignment','Center','VerticalAlignment','top','interpreter', 'latex','FontSize',14);
text(B_tilde,yTicks(1)-0.3,'$\tilde{B}$', 'HorizontalAlignment','Center','VerticalAlignment','top','interpreter', 'latex','FontSize',14); hold on;
text(b_bar,yTicks(1)-0.3,'$\bar{B}$', 'HorizontalAlignment','Center','VerticalAlignment','top','interpreter', 'latex','FontSize',14); hold on;
text(B_exit_num,yTicks(1)-0.3,'$B^{h}$', 'HorizontalAlignment','Center','VerticalAlignment','top','interpreter', 'latex','FontSize',14); hold on;
plot(B_back_l_set,V_d_bound,'b*');
plot(B_back_r_set(1:end-1),V_d_bound(2:end),'g*');
for tt=1:T+1
    grid on; yvals=ylim; line([B_back_set(tt) B_back_set(tt)],[yvals(1) yvals(2)],'Color',[0.1 0.3 0.5],'LineWidth',0.5,'LineStyle',':'); grid minor;
end
for tt=1:N_tau
    grid on; yvals=ylim; line([B_back_l_set(tt) B_back_l_set(tt)],[yvals(1) yvals(2)],'Color',[0.9 0.3 0.5],'LineWidth',0.7,'LineStyle',':'); grid minor;
    grid on; yvals=ylim; line([B_back_r_set(tt) B_back_r_set(tt)],[yvals(1) yvals(2)],'Color',[0.2 0.6 0.5],'LineWidth',0.7,'LineStyle',':'); grid minor;
end
grid on; yvals=ylim; line([B_back_lim B_back_lim],[yvals(1) yvals(2)],'Color',[0.5 0.6 0.6],'LineWidth',1,'LineStyle','-.'); grid minor;


%end
% for tt=1:T+1
%     grid on; yvals=ylim; line([B_back_set(tt) B_back_set(tt)],[yvals(1) yvals(2)],'Color',[0.1 0.3 0.5],'LineWidth',0.5,'LineStyle',':'); grid minor;
% end
% for tt=1:N_tau
%     grid on; yvals=ylim; line([B_back_l_set(tt) B_back_l_set(tt)],[yvals(1) yvals(2)],'Color',[0.9 0.3 0.5],'LineWidth',1,'LineStyle',':'); grid minor;
%     grid on; yvals=ylim; line([B_back_r_set(tt) B_back_r_set(tt)],[yvals(1) yvals(2)],'Color',[0.2 0.6 0.5],'LineWidth',1,'LineStyle',':'); grid minor;
% end
grid on; yvals=ylim; line([B_back_lim B_back_lim],[yvals(1) yvals(2)],'Color',[0.5 0.6 0.6],'LineWidth',1,'LineStyle','-.'); grid minor;
if printit==1
    orient landscape;
    % saveas(gcf,'F_valuefunction','pdf');
    ax= gca;
    exportgraphics(ax,'F_valuefunction_full.pdf','BackgroundColor','none');
end

figure('Name','Backward Sets (Debt Policy)','NumberTitle','off')
% Colors
color_t1=[0.7 0.7 0.7];
color_t0=[0.0 0.0 0.6];
%fplot(@(B) val_high(B),range,'LineWidth',1.5,'LineStyle','-','Color',color_t0); hold on; axis tight;
% fplot(@(B) val_low(B),range,'LineWidth',1.5,'LineStyle','-','Color',color_t1); grid on;
% fplot(@(B) val_high(B)+log(1/q)/(1-beta),range,'LineWidth',3,'LineStyle',':','Color',color_t1); hold on; axis tight;
plot(B_vec,Bp_out,'LineWidth',2,'LineStyle','-','Color','k'); grid on; drawnow; hold on;
plot(B_vec,B_vec,'LineWidth',2,'LineStyle',':','Color','k');
legend('$B''$','$B$','Interpreter','Latex','AutoUpdate','off','Box','Off','Location','northwest')
yvals=ylim; xvals=xlim;
set(gca,'XTick',Bstar,'XTickLabel',[]);
set(gca,'YTick',[yvals(1) yvals(2)],'YTickLabel',[]);
yTicks = get(gca,'ytick');
xTicks = get(gca, 'xtick');
line([xvals(1) xvals(2)],[Bstar Bstar],'Color','k','LineWidth',1,'LineStyle',':');
line([xvals(1) xvals(2)],[B_tilde B_tilde],'Color','k','LineWidth',1,'LineStyle',':');
yvals=ylim; line([Bstar Bstar],[yvals(1) yvals(2)],'Color','k','LineWidth',1,'LineStyle','--');
yvals=ylim; line([Bast Bast],[yvals(1) yvals(2)],'Color','k','LineWidth',1,'LineStyle','--');
yvals=ylim; line([B_tilde B_tilde],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.5],'LineWidth',1,'LineStyle','-.'); % grid minor;
yvals=ylim; line([B_exit_num B_exit_num],[yvals(1) yvals(2)],'Color',[0.5 0.0 0.0],'LineWidth',1,'LineStyle','-.'); % grid minor;
% grid on; yvals=ylim; line([Bhat(B_tilde,B_tilde) Bhat(B_tilde,B_tilde)],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.5],'LineWidth',1,'LineStyle',':');
text(xTicks(1),yTicks(1)-0.3,'$B^{\star}$', 'HorizontalAlignment','Center','VerticalAlignment','top','interpreter', 'latex','FontSize',14);
text(Bast,yTicks(1)-0.3,'$B^{\ast}(\tilde{B})$', 'HorizontalAlignment','Center','VerticalAlignment','top','interpreter', 'latex','FontSize',14);
text(B_tilde,yTicks(1)-0.3,'$\tilde{B}$', 'HorizontalAlignment','Center','VerticalAlignment','top','interpreter', 'latex','FontSize',14); hold on;
text(b_bar,yTicks(1)-0.3,'$\bar{B}$', 'HorizontalAlignment','Center','VerticalAlignment','top','interpreter', 'latex','FontSize',14); hold on;
text(B_exit_num,yTicks(1)-0.3,'$B^{h}$', 'HorizontalAlignment','Center','VerticalAlignment','top','interpreter', 'latex','FontSize',14); hold on;
for tt=1:T+1
    grid on; yvals=ylim; line([B_back_set(tt) B_back_set(tt)],[yvals(1) yvals(2)],'Color',[0.1 0.3 0.5],'LineWidth',0.5,'LineStyle',':'); grid minor;
end
for tt=1:N_tau
    grid on; yvals=ylim; line([B_back_l_set(tt) B_back_l_set(tt)],[yvals(1) yvals(2)],'Color',[0.9 0.3 0.5],'LineWidth',0.7,'LineStyle',':'); grid minor;
    grid on; yvals=ylim; line([B_back_r_set(tt) B_back_r_set(tt)],[yvals(1) yvals(2)],'Color',[0.2 0.6 0.5],'LineWidth',0.7,'LineStyle',':'); grid minor;
end
grid on; yvals=ylim; line([B_back_lim B_back_lim],[yvals(1) yvals(2)],'Color',[0.5 0.6 0.6],'LineWidth',1,'LineStyle','-.'); grid minor;
if printit==1    
    orient landscape;
  %  set(gca,'LooseInset',get(gca,'TightInset'));
  %  saveas(gcf,'F_bprime','pdf');
    ax = gca;
    exportgraphics(ax,'F_bprime_full.pdf','BackgroundColor','none');
end

figure('Name','Backward Sets (C Policy)','NumberTitle','off')
% Colors
color_t1=[0.7 0.7 0.7];
color_t0=[0.0 0.0 0.6];
%fplot(@(B) val_high(B),range,'LineWidth',1.5,'LineStyle','-','Color',color_t0); hold on; axis tight;
% fplot(@(B) val_low(B),range,'LineWidth',1.5,'LineStyle','-','Color',color_t1); grid on;
% fplot(@(B) val_high(B)+log(1/q)/(1-beta),range,'LineWidth',3,'LineStyle',':','Color',color_t1); hold on; axis tight;
% plot(B_vec,C_high(B_vec),'LineWidth',2,'LineStyle',':','Color','k'); 
plot(B_vec,C_out,'LineWidth',2,'LineStyle','-','Color','k'); drawnow; grid on; hold on;
legend('$C(B;\tilde{B}>0)$','Interpreter','Latex','AutoUpdate','off','Box','Off')
yvals=ylim; xvals=xlim;
set(gca,'XTick',Bstar,'XTickLabel',[]);
set(gca,'YTick',[yvals(1) yvals(2)],'YTickLabel',[]);
yTicks = get(gca,'ytick');
xTicks = get(gca, 'xtick');
grid on; yvals=ylim; line([Bstar Bstar],[yvals(1) yvals(2)],'Color','k','LineWidth',1,'LineStyle','--');
%grid on; yvals=ylim; line([Bast Bast],[yvals(1) yvals(2)],'Color','k','LineWidth',1,'LineStyle','--');
grid on; yvals=ylim; line([B_tilde B_tilde],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.5],'LineWidth',1,'LineStyle','-.'); grid minor;
% grid on; yvals=ylim; line([B_exit B_exit],[yvals(1) yvals(2)],'Color',[0.5 0.0 0.0],'LineWidth',1,'LineStyle','-.'); grid minor;
grid on; yvals=ylim; line([B_exit_num B_exit_num],[yvals(1) yvals(2)],'Color',[0.5 0.0 0.0],'LineWidth',1,'LineStyle','-.'); grid minor;
% for TT=(T_trans:-1:2)
%     scatter(B_path(TT),c_path(TT),'r*');
% end
% grid on; yvals=ylim; line([Bhat(B_tilde,B_tilde) Bhat(B_tilde,B_tilde)],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.5],'LineWidth',1,'LineStyle',':');
text(xTicks(1),yTicks(1)-0.1,'$B^{\star}$', 'HorizontalAlignment','Center','VerticalAlignment','bottom','interpreter', 'latex');
%text(Bast,yTicks(1)-0.1,'$B^{\ast}(\tilde{B})$', 'HorizontalAlignment','Center','VerticalAlignment','bottom','interpreter', 'latex');
text(B_tilde,yTicks(1)-0.1,'$\tilde{B}$', 'HorizontalAlignment','Center','VerticalAlignment','bottom','interpreter', 'latex'); hold on;
for tt=1:T+1
    grid on; yvals=ylim; line([B_back_set(tt) B_back_set(tt)],[yvals(1) yvals(2)],'Color',[0.1 0.3 0.5],'LineWidth',0.5,'LineStyle',':'); grid minor;
end
for tt=1:N_tau
    grid on; yvals=ylim; line([B_back_l_set(tt) B_back_l_set(tt)],[yvals(1) yvals(2)],'Color',[0.9 0.3 0.5],'LineWidth',0.7,'LineStyle',':'); grid minor;
    grid on; yvals=ylim; line([B_back_r_set(tt) B_back_r_set(tt)],[yvals(1) yvals(2)],'Color',[0.2 0.6 0.5],'LineWidth',0.7,'LineStyle',':'); grid minor;
end
grid on; yvals=ylim; line([B_back_lim B_back_lim],[yvals(1) yvals(2)],'Color',[0.5 0.6 0.6],'LineWidth',1,'LineStyle','-.'); grid minor;
if printit==1    
    orient landscape;
  %  set(gca,'LooseInset',get(gca,'TightInset'));
  %  saveas(gcf,'F_bprime','pdf');
    ax = gca;
    exportgraphics(ax,'F_cons_full.pdf','BackgroundColor','none');
end

%% Main Simulations in Partial Equilibrium
%% Simulated Paths
% Simulation settings
T = 25;  % horizon

% Initial index positions (grid-aligned!)
init_bs = [5, N_b-10, N_b - 10, N_b - 150];

colors = lines(length(init_ks));
B_hist=[]; K_hist=[]; C_hist=[];
for i = 1:length(init_ks)
    k_idx = init_ks(i);
    b_idx = init_bs(i);

    K_path = zeros(T, 1);
    B_path = zeros(T, 1);
    C_path = zeros(T, 1);
    for t = 1:T
        K_path(t) = K_vec(k_idx);
        B_path(t) = B_vec(b_idx);
       
        % Advance to next state using grid indices
        K_next = Kp_out(k_idx, b_idx);
        B_next = Bp_out(k_idx, b_idx);

        % Update Consumption
        C_path(t) = C_out(k_idx, b_idx);

        % Find indices for next step
        [~, k_idx] = min(abs(K_vec - K_next));
        [~, b_idx] = min(abs(B_vec - B_next));
    end
    plot(B_path, K_path, '-', 'Color', colors(i,:), 'LineWidth', 1.5);
    plot(B_path(end), K_path(end), 'o', 'MarkerSize', 3, ...
         'Color', colors(i,:), 'MarkerFaceColor', colors(i,:));
        plot(B_path, K_path, '<', 'MarkerSize', 3, ...
         'Color', colors(i,:), 'MarkerFaceColor', colors(i,:));

   % Update histories
   B_hist=[B_hist B_path];
   K_hist=[K_hist K_path];
   C_hist=[C_hist C_path];
end
axis tight;
legend(arrayfun(@(i) sprintf('Init %d', i), 1:length(init_ks), 'UniformOutput', false), 'Location','best');
line([B_tilde B_tilde],[K_vec(1) K_vec(end)],'LineWidth',2,'Color','r');
line([Bstar Bstar],[K_vec(1) K_vec(end)],'LineWidth',2,'Color','r');

figure('Name','B_simulation');
plot(1:T,B_hist)
xlabel('Time'); ylabel('Debt');
grid on; axis tight;

figure('Name','K_simulation');
plot(1:T,K_hist)
xlabel('Time'); ylabel('Capital');
grid on; axis tight;

figure('Name','C_simulation');
plot(1:T,C_hist)
xlabel('Time'); ylabel('Consumption');
grid on; axis tight;
