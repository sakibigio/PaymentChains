%% Special Solutions
% Solution without payment networks: exogenous prices
funcs.val_high=@(B) log((1-(1-params.beta)*B))/(1-params.beta)  ; % no chained consumption
funcs.val_low=@(B) log((1-(1-params.beta)*B)/params.q)/(1-params.beta) ; % all consumption is at chained price q
funcs.C_low=@(B) (1-(1-params.beta)*B)/params.q; % annuity value of consumption

%% Value Function Iteration
% define functions - B_p is future debt, determining expenditures today
funcs.E     = @(B_p,B) 1-B+B_p/params.R                    ;

% knowing expenditures and credit line, you can determine optimal expenditures in spot and chained, and determine total consumption
funcs.S_w   = @(E,B,B_tilde) min(max(B_tilde-B,0),E);
funcs.X_w   = @(E,B,B_tilde) (E-funcs.S_w(E,B,B_tilde))/params.q ;
funcs.C_w   = @(E,B,B_tilde) funcs.X_w(E,B,B_tilde)+funcs.S_w(E,B,B_tilde)          ;
 
%% Start Figures
% figure
% fplot(@(B) val_high(B),range,'LineWidth',3,'LineStyle','-'); hold on; axis tight;
% fplot(@(B) val_low(B),range,'LineWidth',3,'LineStyle','-'); grid on;% fplot(@(B) val_high(B)+log(1/q)/(1-beta),range,'LineWidth',3,'LineStyle',':'); hold on; axis tight;
