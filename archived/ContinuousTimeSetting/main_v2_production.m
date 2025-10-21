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
plotiter=1;

%% Optimization Conditions (create subroutine)
tol_dyn = 1e-8;
tol     = 1e-4;
cond    = 2*tol;
options = optimset('TolFun',1e-6,'Display','iter'); % At Optimization
%options_dyn=optimset('Algorithm',{'levenberg-marquardt',0.01},'MaxFunEvals',10*T,'MaxIter',400,'TolFun',tol_dyn,'Display','iter');
% options_dyn=optimset('Algorithm','trust-region-dogleg','MaxFunEvals',1000,'MaxIter',400,'TolFun',tol_dyn,'Display','iter');

%% Parameters (create subroutine)
% Model Parameters - Preferences and Technology
freq    = 1     ; % 1 year over freqency
gamma   = 1     ; % agent's risk aversion
rho     = 0.0281; % agent's individual discount rates; target: 1% (quarterly)
q       = 2     ; % high-price values

% Bliss point
c_bliss = 0; 

% Production Function
alpha=0.5; % DRS
f=@(k) k.^(alpha)  ;
f_inv=@(x) x.^(1/alpha);
f_p=@(k) alpha*k.^(alpha-1); % f prime
f_pinv=@(f) (f/alpha).^(1/(alpha-1)); % f prime
pi_opt=@(q) f(f_pinv(q))-q.*f_pinv(q);

% Borrowing Limit - nominal
pi_1=pi_opt(q);

scale=2; % -> bad coding! move up
s_bl  = -1/scale*0.99*pi_1/rho       ;  % "soft constraint"
s_max  = (pi_1/rho)*1  ;  % debt limit
s_bar = scale*s_bl;

% Vectorizations
mu_w    = (pi_1)';
% sigma_w = [s1 s2]';

% [1.ba] Guesses
% rs_o  = 0.01;%rho-rsp_ss-2e-2;
% T_o   = 1/2;

% [IV] Functional Forms 
% Utility Functions
if gamma~=1
    URF=@(c) ((c-c_bliss).^(1-gamma)-1)/(1-gamma); % Utility Return Function
    Uprime=@(c) (c-c_bliss).^(-gamma); % Utility Return Function
else
    URF=@(c)   log(c-c_bliss); % Utility Return Function
    Uprime=@(c) (c-c_bliss).^(-1); % Utility Return Function
end
Upinv=@(dv) dv.^(-1/gamma)+c_bliss;              % Utility Return Function

% Saving parameters
parameters.gamma=gamma;



%% Discretization Preferences
% Approximate Amount of grid points
N     = 2000      ;  % Number of Gridpoints in Real Wealth Space

% State-Space Grid
dt     = 1/12     ; % Double check code for dt.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Discretization Block
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Savings space and consumption space Gus
s_vec      = linspace(s_bar, s_max,N)'        ; 
ds         = s_vec(2)-s_vec(1); ds_check=(s_max-s_bar)/(N-1);
s_bl_index = (s_vec<=s_bl);

% Compute Savers
index_o=find(s_vec>=0,1,'first')          ;

% Initialized Variables
% mu_ss    = zeros(1,N); % Drift at Steady State
% sigma_ss = zeros(1,N); % Volatility at Steady State

% addpath(genpath([cd '/Codes2']));

% Variables defined as global come in
% MPCC_globals;

% Initial Guess for consumption vector
% p  = inputs.p;
rs      = rho;

% Ttransf = X(2);
% Ttransf=0;
% transratio = X(2);

% [1] Solve for Steady State for each case
% Construct Interest Rate Vector
r_vec=ones(N,1)*rs  ; % Initialize i
% r_vec(s_vec<=0)=rs; % pin down value
q_vec=ones(N,1); 
q_vec(s_bl_index)=q; 
w1=pi_opt(q_vec);

% Initial conjecture for value
c_guess = (w1+r_vec.*s_vec)./q_vec;
distance= s_vec-min(s_vec);
weights = exp(-0.1*distance)./(0.5+0.5*exp(-0.1*distance));
V       = ((1-weights).*URF(w1+r_vec.*s_vec)+(weights).*URF((w1+r_vec.*s_vec)/q))/rho;
V_bar   = URF(max(w1)+r_vec.*s_vec)/rho;
V_ubar  = URF((min(w1)+r_vec.*s_vec)/q)/rho;

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
    cond=max(abs(V./V_in-1));
    cond

    % Plot if needed
    if plotiter==1
        plot(s_vec,V); drawnow; hold on;
    end
    iter=iter+1; 
end

% Plot Value Function
figure('Name','Value Function')
plot(s_vec,V); grid on; hold on;
plot(s_vec,V_bar);
plot(s_vec,V_ubar);
 
% Colecting Steady State Objects
c_ss     = c     ;
muF_ss   = muF   ;
muB_ss   = muB   ;
V_ss     = V     ;



%% Analytic Solution
% % Solve v_bar
% c_star=w1+rs*s_bl;
% RHS_vbar=@(v) (rho+log(q)+log(c_star)+log(v))/(w1+rho*s_bl);
% fplot(@(x) RHS_vbar(x),[0,2]);v_bar=fsolve(@(v) RHS_vbar(v)-v,1);
% 
% % Compare these with paper and pencil solution. How do they look?
% % solve particular case:
% 
% c_back=fsolve(@(c) Uprime(c)/q-(URF(c_star)-URF(c))/(w1-q*c+rho*s_bl),w1/25);
% c_drop=c_back;
% 
% figure('Name','Consumption')
% plot(s_vec,c); hold on;
% plot(s_vec,c*0+c_back,'k--'); hold on;
% axis tight; grid on;







%% Steady State Solution
function res=MPCC_solve_steady(X,inputs)
% Variables defined as global come in
MPCC_globals;

% Initial Guess for consumption vector
p  = inputs.p;
rs      = X(1);
% Ttransf = X(2);
% Ttransf=0;


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