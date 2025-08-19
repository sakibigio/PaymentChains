% Continuout-Time Solutions - Parameters
freq    = 1     ; % 1 year over freqency
gamma   = 1     ; % agent's risk aversion
rho     = 0.051; % agent's individual discount rates; target: 1% (quarterly)
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