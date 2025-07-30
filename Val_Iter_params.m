%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%
% Value functions - Partial Equilibrium
params.beta= 0.8   ; % 0.8
params.R   = 1/params.beta;
params.q   = 1.75  ; % 1.75
params.b_bar=1/(1-params.beta); % one unit of labor input, or annuity value: natural borrowing limit
params.B_tilde=0.4*params.b_bar; % 0.4
params.Bstar=1/params.beta*(params.B_tilde-1); % Threshold
params.Bast  = params.B_tilde/(params.q^(-1)-(params.q^(-1)-1)*(1-params.beta)*params.B_tilde);


