% Main Functional forms
C_s   = @(B) (1-params.beta).*B;
S_w   = @(B,B_tilde) min(max(B_tilde-B,0),1-(1-params.beta)*B);
mu_Z  = @(B,B_tilde) 1-C_s(B)-S_w(B,B_tilde);
q_Z   = @(B,B_tilde) q_mu(mu_Z(B,B_tilde),params.delta);
X_w   = @(B,B_tilde) mu_Z(B,B_tilde)./q_Z(B,B_tilde);
Y_Z   = @(B,B_tilde) 1-mu_Z(B,B_tilde)+mu_Z(B,B_tilde)./q_mu(mu_Z(B,B_tilde),params.delta);

% Derived Functions
% Average price
Q_Z   = @(B,B_tilde) (1-(1-params.beta)*B)./(X_w(B,B_tilde)+S_w(B,B_tilde));

% Marginal Price of Savings
q_B=@(B,B_tilde) (1+(q_Z(B,B_tilde)-1).*(B_tilde<B));

% Marginal Price of eExpenditures
q_E=@(B,B_tilde) (1+(q_Z(B,B_tilde)-1).*(B_tilde<1+params.beta*B));

% Marginal Inflation
Pi =@(B,B_p,B_tilde,B_tilde_p) q_B(B_p,B_tilde_p)./q_E(B,B_tilde);

% Ratio of Average to Marginal Price (expenditure and savings)
Qq_E=@(B,B_tilde) Q_Z(B,B_tilde)./q_E(B,B_tilde);
Qq_B=@(B,B_tilde) Q_Z(B,B_tilde)./q_B(B,B_tilde);

% Euler_rhs=@(B,B_p,B_tilde,B_tilde_p)  B .* (1-params.beta*B_p)./(1-params.beta*B).*Q_Z(B,B_tilde)./Q_Z(B_p,B_tilde_p).*Pi(B,B_p,B_tilde,B_tilde_p);
% Euler_rhs=@(B,B_p,B_tilde,B_tilde_p)  B./(1-(1-params.beta)*B).*Q_Z(B,B_tilde).*(X_w(B_p,B_tilde_p)+S_w(B_p,B_tilde_p)).*Pi(B,B_p,B_tilde,B_tilde_p);
Euler_t_p = @(B_p,B_tilde_p) 1./(1/B_p-(1-params.beta)).*Qq_B(B_p,B_tilde_p);
Euler_t   = @(B,B_tilde) 1./(1/B-(1-params.beta)).*Qq_E(B,B_tilde);

B_p_res=@(B,B_p,B_tilde,B_tilde_p) Euler_t_p(B_p,B_tilde_p)  - Euler_t(B,B_tilde);

%B_p_res=@(B,B_p,B_tilde,B_tilde_p) (1/B-(1-params.beta))./Qq_E(B,B_tilde)-(1/B_tilde-(1-params.beta))./Qq_B(B_tilde,B_tilde_p);

% Critical Points
Bstar = @(B_tilde) 1/params.beta*(B_tilde-1);
Bast  = @(B_tilde_p) B_tilde_p/(q_Z(B_tilde_p,B_tilde_p).^(-1)-(q_Z(B_tilde_p,B_tilde_p).^(-1)-1)*(1-params.beta)*B_tilde_p);
Bhat  = @(B_tilde,B_tilde_p) fsolve(@(x) Euler_t(x,B_tilde)-1./(1/Bstar(B_tilde)-(1-params.beta)),Bstar(B_tilde)+0.2); 

%% Key solutions for planner
P    = @(B,B_tilde,theta) (1-theta)*log((1-params.beta)*B)+theta*log(S_w(B,B_tilde)+X_w(B,B_tilde));
B_Ram= @(B_tilde,theta) fminsearch(@(x) -P(x,B_tilde,theta),0.01);

%% Checks: all satisfied
% make sure C+S+X adds to Y
y_check=@(B,B_tilde) Y_Z(B,B_tilde)-(C_s(B)+S_w(B,B_tilde)+X_w(B,B_tilde));
income_check=@(B,B_tilde) 1-(C_s(B)+S_w(B,B_tilde)+q_Z(B,B_tilde).*X_w(B,B_tilde));

% Saving in structure
funcs.y_check=y_check;
funcs.income_check=income_check;
funcs.C_s=C_s;
funcs.S_w=S_w;
funcs.mu_Z=mu_Z;
funcs.q_Z=q_Z;
funcs.X_w=X_w;
funcs.Y_Z=Y_Z;
funcs.Q_Z=Q_Z;
funcs.Qq_E=Qq_E;
funcs.Qq_B=Qq_B;
funcs.Bstar=Bstar;
funcs.Bast=Bast;
funcs.B_p_res=B_p_res;
funcs.B_Ram=B_Ram;
funcs.P=P;