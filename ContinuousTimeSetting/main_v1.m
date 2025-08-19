%% Main Code for Continuous-Time Version of Payment-Chain Model
% Uses Implicit Finite Difference
% (c) Bigio et al
%
% October 2025
% Main block happens here
% written by S. Bigio
clear; close all;
%  nameplot = 'rsp';

%% Running Preferences
plotiter=0;

%% Optimization Conditions (create subroutine)
tol_dyn = 1e-7;
tol     = 1e-10;
options = optimset('TolFun',1e-6,'Display','iter'); % At Optimization
%options_dyn=optimset('Algorithm',{'levenberg-marquardt',0.01},'MaxFunEvals',10*T,'MaxIter',400,'TolFun',tol_dyn,'Display','iter');
% options_dyn=optimset('Algorithm','trust-region-dogleg','MaxFunEvals',1000,'MaxIter',400,'TolFun',tol_dyn,'Display','iter');

%% Parameters (create subroutine)
% Model Parameters - Preferences and Technology
freq    = 1     ; % 1 year over freqency
gamma   = 1     ; % agent's risk aversion
rho     = 0.0351; % agent's individual discount rates; target: 1% (quarterly)
w1      = 1     ; % income

% Prices 
rs      = rho; % risk free rate
q       = 2  ; % Exogenous price

% Bliss point
c_bliss = 0; 

% Borrowing Limit - nominal
scale=2; % -> bad coding! move up
s_bl  = -1/scale*0.99*w1/rho       ;  % "soft constraint"
s_max  = (w1/rho)*0.1  ;  % debt limit
s_bar = scale*s_bl;

% Vectorizations
mu_w    = (w1)';

% Functional Forms 
% Utility Functions
if gamma~=1
    URF=@(c) ((c-c_bliss).^(1-gamma)-1)/(1-gamma); % Utility Return Function
    Uprime=@(c) (c-c_bliss).^(-gamma); % Utility Return Function
    Uinv=@(v) ((1-gamma)*v).^(1/1-gamma)+c_bliss;              % Utility Return Function
else
    URF=@(c)   log(c-c_bliss); % Utility Return Function
    Uprime=@(c) (c-c_bliss).^(-1); % Utility Return Function
    Uinv=@(v) exp(v)+c_bliss;              % Utility Return Function
end
Upinv=@(dv) dv.^(-1/gamma)+c_bliss;              % Utility Return Function


% Saving parameters
parameters.gamma=gamma;

%% Discretization Preferences
% Approximate Amount of grid points
N     = 5000      ;  % Number of Gridpoints in Real Wealth Space

% State-Space Grid
dt     = 1/60     ; % Double check code for dt.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Discretization Block
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Savings space and consumption space Gus
s_vec      = linspace(s_bar, s_max,N)'        ; 
ds         = s_vec(2)-s_vec(1); ds_check=(s_max-s_bar)/(N-1);
s_bl_index = (s_vec<s_bl);
s_unc_index= (s_vec>=s_bl);

% Compute Savers
index_o=find(s_vec>=0,1,'first')          ;

% [1] Solve for Steady State for each case
% Construct Interest Rate Vector
r_vec=ones(N,1)*rs  ; % Initialize i
q_vec=ones(N,1)     ; 
q_vec(s_bl_index)=q ; 

% Initial conjecture for value
c_guess = (w1+r_vec.*s_vec)./q_vec;
distance= s_vec-min(s_vec);
weights = exp(-0.1*distance)./(0.5+0.5*exp(-0.1*distance));
V       = ((1-weights).*URF(w1+r_vec.*s_vec)+(weights).*URF((w1+r_vec.*s_vec)/q))/rho;
V_bar   = URF(w1+r_vec.*s_vec)/rho;
V_ubar  = URF((w1+r_vec.*s_vec)/q)/rho;

figure('Name','Guess')
plot(s_vec,V); grid on; hold on;
plot(s_vec,V_bar);
plot(s_vec,V_ubar);

%% Main Step For solution
cond=2*tol;
iter = 0;
dvF = zeros(N,1); dvB = zeros(N,1);
while cond>tol % value function iteration
    V_in = V;
    
    % SB note: document here:
    ss   = find(s_bl_index,1,'last');
    HJBupdate;    
    
    [V,A]=HJB_implicit(V_in,U,s_vec,muF,muB,rho,Ib,If,dt);    
    
    % Update condition
    % cond=max(abs(V./V_in-1));
    cond=max(abs(Uinv(rho*V)./Uinv(rho*V_in)-1));
    % display(['Error condition: ' num2str(cond)]);

    % Plot if needed
    if plotiter==1
        plot(s_vec,V); drawnow; hold on;
    end
    iter=iter+1; 
end

% Plot Value Function
figure('Name','Value Function')
plot(s_vec,V_bar,'LineWidth',1);  grid on; hold on;
plot(s_vec,V_ubar,'LineWidth',1);
plot(s_vec,V,'LineWidth',2,'LineStyle','-.');

% Colecting Steady State Objects
c_ss     = c     ;
muF_ss   = muF   ;
muB_ss   = muB   ;
V_ss     = V     ;

%% Analytic Solution
% solve consumption at deleveraging phase:
c_star=w1+rs*s_bl;
c_back=fsolve(@(c) Uprime(c)/q-(URF(c_star)-URF(c))/(w1-q*c+rho*s_bl),w1/25);
c_drop=c_back;
v_bar = (URF(c_star)-URF(c_back))/(w1-q*c_back+rho*s_bl);

% Solve v_bar
RHS_b_h=@(b_h) URF((w1+rho*b_h)/q);
LHS_b_h=@(b_h) URF(c_back)+v_bar*(w1-q*c_back+rho*b_h);


b_h=fsolve(@(x) RHS_b_h(x)-LHS_b_h(x),s_bl);
index_aux=((s_vec>b_h)&(s_vec<s_bl));
index_aux2=((s_vec<=b_h));
plot(s_vec(index_aux),URF((w1+rho*b_h)/q)/rho+v_bar*(s_vec(index_aux)-b_h),'LineWidth',1,'LineStyle',':','Color','r');
plot(s_vec(s_bl_index),URF((w1+rho*s_vec(s_bl_index))/q)/rho,'LineWidth',1,'LineStyle',':','Color','r');
plot(s_vec(s_unc_index),URF((w1+rho*s_vec(s_unc_index)))/rho,'LineWidth',1,'LineStyle',':','Color','r');

figure
fplot(@(x) RHS_b_h(x),[-30,30]); hold on
fplot(@(x) LHS_b_h(x),[-30,30]); hold on

figure('Name','Consumption')
plot(s_vec(index_aux),c(index_aux),'Color','b'); hold on; grid on;
plot(s_vec(index_aux2),c(index_aux2),'Color','b'); hold on; grid on;
plot(s_vec(s_unc_index),c(s_unc_index),'Color','b'); 
line([s_vec(1) s_vec(end)],[c_back c_back],'LineStyle',':','LineWidth',1); 
scatter(s_bl,c_star,40,'MarkerFaceColor','b')
scatter(s_bl,c_back,40,'MarkerFaceColor','w')
scatter(b_h,(w1+rho*b_h)/q,40,'MarkerFaceColor','b')
scatter(b_h,c_back,40,'MarkerFaceColor','w')
ylimits=ylim;
line([b_h b_h],[ylimits(1) ylimits(2)],'LineWidth',1);
line([s_bl s_bl],[ylimits(1) ylimits(2)],'LineWidth',1);
axis tight;

%% Steady State Solution
function res=MPCC_solve_steady(X,inputs)
% Variables defined as global come in
MPCC_globals;

% Initial Guess for consumption vector
p  = inputs.p;
rs      = X(1);

% [1] Solve for Steady State for each case
% Construct Interest Rate Vector
r_vec=s_vec*0+rs       ; % Initialize i
r_vec(s_vec<=0)=rs+rsp_ss; % pin down value

% Initial conjecture for value
c_guess = w1+r_vec.*s_vec;
V       = URF(c_guess)/rho;

% Calls the Solver
MPCC_steady;

% Pick a clearing condition to solve
if strcmp(clearcond,'Y');
    res=Z_Y_ss;
elseif strcmp(clearcond,'S');
    res=Z_S_ss;
end
% res(2)  = Z_T_ss;

end

%% HJB_implicit
function [V,A]=HJB_implicit(V_in,U,s,muF,muB,rho,Ib,If,Delta)
    
    % Vectorization
    N = length(s)         ; % Length
    ds = s(2)-s(1);
    
    % Plus diagonal
    MU  = If.*muF;
    DU  = (MU)/ds   ; % Negative from right of HJB 
    
    % Minus Diagonal Terms
    MD = -Ib.*muB;
    DD = (MD)/ds     ; % Negative from right of HJB
    
    % On Diagonal Terms
    D0 = -DU-DD;
    
    % Making it Compatible with Sparse Matrix Notation
    DU = [0; DU(1:end-1)];
    DD = [DD(2:end); 0];
    
    A  = spdiags([DD D0 DU],-1:1,N,N);
    
    V = (speye(N)*(1/Delta+rho)-A)\(U+V_in/Delta);

end