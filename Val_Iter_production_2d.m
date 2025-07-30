%% Payment Chains
% (c) 
% tasks:
% write it in terms of assets, not debt
% write it with i, do not obtain all values but optimize for speed
% make sure distortion in investment is real
 
% version with only consumption, no production. 
clear; close all;
tol=10e-6;
printinfo=1;
printit=0;
plotit=1;

%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%
% Grid Points
N_b=60; % grid points for debt level
N_k=60; % grid points for capital

% Value functions - Partial Equilibrium
beta= 0.8           ; % 0.8
R   = 1/beta        ;
q   = 1.75          ; % 1.75

% Production Function
alpha=0.5; % DRS
f=@(k) k.^(alpha)  ;
f_inv=@(x) x.^(1/alpha);
f_p=@(k) alpha*k.^(alpha-1); % f prime
f_pinv=@(f) (f/alpha).^(1/(alpha-1)); % f prime

% We discussed two possible models, one where agents buy capital in
% the market and one where they don't
% Version: capital bought in the market:
K_high=f_pinv(R)  ;
K_low =f_pinv(R*q);

% Special Forms
b_bar=(f(K_low)-q*K_low)/(1-beta)    ; % one unit of labor input, or annuity value: natural borrowing limit
B_tilde=0.4*b_bar   ; % 0.4
Bstar=1/beta*(B_tilde-1); % Threshold
Bast  = B_tilde/(q^(-1)-(q^(-1)-1)*(1-beta)*B_tilde);

% Support of values of debt
%a_max=f(K_high)-q*K_high;
%a_bar=f(K_low)-R*b_bar;
%range=[a_bar+(10e-10) a_max*2]    ;
%A_vec=linspace(range(1),range(end),1000); % grid space
K_vec=linspace(K_low*0.5,K_high,N_k);
B_vec=linspace(0,b_bar-1e-6,N_b);

%% Functional Forms
% Extreme values - Initial Guesses - Capital at Steady State pay annuity:
C_high=@(K,B)   (f(K)-K)*ones(1,length(B))-ones(length(K),1)*(B-B/R); % annuity value of consumption
val_high=@(K,B) log(C_high(K,B))/(1-beta)  ; % no chained consumption
C_low=@(K,B)    ((f(K)-q*K)*ones(1,length(B))-ones(length(K),1)*(B-B/R)); % annuity value of consumption
val_low=@(K,B)  log(C_low(K,B))/(1-beta) ; % all consumption is at chained price q

rango_plot=[K_low K_high 0 b_bar];
if plotit==1
    figure('Name',"Value Function Bounds")
    fsurf(@(K,B) val_high(K,B),rango_plot); hold on; axis tight;
    fsurf(@(K,B) val_low(K,B),rango_plot); grid on;% 
    figure('Name',"Consumption")
    fsurf(@(K,B) C_high(K,B),rango_plot); hold on; axis tight;
    fsurf(@(K,B) C_low(K,B),rango_plot); grid on;% 
end

%% Value Function Iteration
% define functions - B_p stands for B_prime (future debt), determining expenditures today
E     = @(B_p,K,B) f(K)-B+B_p/R                    ; % total expenditures as functions of current wealth minus debt

% knowing expenditures and credit line, you can determine optimal expenditures in spot and chained, and determine total consumption
S_w   = @(E,B,B_tilde) min([max(B_tilde-B,0)*ones(1,length(E)); E]);
X_w   = @(E,B,B_tilde) (E-S_w(E,B,B_tilde))/q ;
C_w   = @(E,B,B_tilde,K_p) ones(length(K_p),1)*(X_w(E,B,B_tilde)+S_w(E,B,B_tilde))-K_p*ones(1,length(E)); % Consumption is goods purchased minus i investment

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Iteration 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize values
V0=log(C_low(K_vec'*0+K_low,B_vec))/(1-beta); 
V_out=V0; C_out=0*V_out; S_out=0*V_out; X_out=0*V_out; Bp_out=0*V_out; Kp_out=Bp_out; i_out=0*V_out;
v_out_ingrid=0*V_out; bp_ingrid=0*V_out;
val_gap=2*tol; % reset tolerance

% Contraction Mapping
figure
while val_gap>tol
    if printinfo
        surf(B_vec,K_vec',V0,'FaceAlpha',0.5); drawnow; hold on;
        disp(['Val at iter:' num2str(val_gap)]);
    end

    % Begin loop among all grid values
    kk=0;
    for K=K_vec
        kk=kk+1;

        % Update definition of RHS of Value Functions
        bb=0;
        for B=B_vec
            bb=bb+1;
            
                % Build grid of values:
                Tv=@(K_p,B_p) log(C_w(E(B_p,K,B),B,B_tilde,K_p'))+beta*V0;
    
                % Solution
                Tv_mat=Tv(K_vec,B_vec);

                % Logical mask for real entries
                realMask = imag(Tv_mat) == 0 & isfinite(Tv_mat);
                
                % Extract real entries
                realVals = Tv_mat(realMask);
                
                % No real values means empty set
                if isempty(realVals)
                    error('No real-valued entries found in Tv_mat.');
                end
                
                % Find the maximum among real values
                [maxRealVal, realLinearIdx] = max(realVals);
                
                % Map back to original matrix indices
                % Find all real entry indices
                realLinearIndices = find(realMask);
                originalLinearIdx = realLinearIndices(realLinearIdx);
                [rowIdx, colIdx] = ind2sub(size(Tv_mat), originalLinearIdx);
                
                % evaluate at optimum to obtain solution
                V_out(kk,bb)=maxRealVal;
                
                % update solutions
                Bp_out(kk,bb)= B_vec(colIdx);
                Kp_out(kk,bb)= K_vec(rowIdx);
                C_out(kk,bb) = C_w(E(Bp_out(kk,bb),K,B),B,B_tilde,Kp_out(kk,bb));
                S_out(kk,bb) = S_w(E(Bp_out(kk,bb),K,B),B,B_tilde);
                X_out(kk,bb) = X_w(E(Bp_out(kk,bb),K,B),B,B_tilde);
        end
    end
    val_gap=abs(max(max(V_out-V0)))/abs(max(max(V0)));
    V0=V_out;
end

%% Plotting Steady-State Solutions:
figure('Name','Value Function')
surf(B_vec,K_vec',V_out,'FaceAlpha',0.5); drawnow; hold on;

figure('Name','Debt Prime')
surf(B_vec,K_vec',Bp_out,'FaceAlpha',0.5); drawnow; hold on;

figure('Name','Capital Prime')
surf(B_vec,K_vec',Kp_out,'FaceAlpha',0.5); drawnow; hold on;

figure;
quiver(B_vec,K_vec,Bp_out,Kp_out, 'k');