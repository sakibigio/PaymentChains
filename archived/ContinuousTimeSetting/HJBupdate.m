

% pre-allocate
dvF = zeros(N,1); dvB = zeros(N,1);

% [1] Update HJB equation
% "steady" value functions U(annuity)/rho. Thus, 
% value function derivative is U'(annuity)*r/rho/q.
% Finite difference approximation -> careful here
dvF(1:end-1) = (V(2:end) - V(1:end-1))./(s_vec(2:end) - s_vec(1:end-1)); % Forward Derivative 
dvF(N)       = Uprime(w1(end) + r_vec(end)*s_vec(end)-c_bliss)*rs/rho; % at final node:
dvB(2:end)   = (V(2:end) - V(1:end-1))./(s_vec(2:end) - s_vec(1:end-1)); % Backward Derivative
dvB(1)       = Uprime((w1(1) + r_vec(1)*s_vec(1))/q-c_bliss)*rs/rho/q; % use Envelope

dvF = max(real(dvF),1e-3);
dvB = max(real(dvB),1e-3);

% Backout Consumption - it is a function of the location
cF = Upinv(dvF)./q_vec           ; % First-order condition
cB = Upinv(dvB)./q_vec           ; % First-order condition
% cF   = min(cF./q_vec,c_bl./q_vec); % include borrowing constraint limit
% cB   = min(cB./q_vec,c_bl./q_vec); % include borrowing constraint limit    
    
% Update drift and vol
muF = zeros(N,1)   ;      % Initialize Vector of Drifts
muB = zeros(N,1)   ;      % Initialize Vector of Drifts

% Update Drifts and Volatilities  
% Assuming we work with risk-free technology at borrowing constraint 
muF       = (r_vec.*s_vec + (w1 - q_vec.*cF));
% muF(1:ss) = r_vec(1:ss).*s_vec(1:ss) + (y1 - cF(1:ss))*p ;                                                             
HcF       = URF(cF) + dvF.*muF;      % forward Hamiltonian

muB       = (r_vec.*s_vec + (w1 - q_vec.*cB));
% muB(1:ss) = r_vec(1:ss).*s_vec(1:ss) + (y1 - cB(1:ss))*p ;                                                             
HcB       = URF(cB) + dvB.*muB;      % backward Hamiltonian

c0        = (r_vec.*s_vec + w1)./q_vec; 
%c0(1:ss)  = r_vec(1:ss).*s_vec(1:ss)/p + w1;
% c0(1:ss)  = min(c0(1:ss),c_bl(1:ss)); % include borrowing constraint limit
dv0       = Uprime(c0);    
H0        = URF(c0);

% Handling with non-convexities - Viscosity Solution
Ineither = (1-(muF>0)) .* (1-(muB<0));
Iunique  = (muB<0).*(1-(muF>0)) + (1-(muB<0)).*(muF>0);
Iboth    = (muB<0).*(muF>0);
Ib       = Iunique.*(muB<0) + Iboth.*(HcB>=HcF);
If       = Iunique.*(muF>0) + Iboth.*(HcF>=HcB);
I0       = Ineither;    
If(end)  = 0; Ib(end) = 1; I0(end)=0;


c    = cF.*If + cB.*Ib + c0.*I0;
adot = muF.*If + muB.*Ib;

% sigma             = ones(N,1)*s1 ; % Initial Vector of Volatilities
% sigma(ss+1:end-1) = s2*p         ; % volatility is real  
        
U =URF(c);