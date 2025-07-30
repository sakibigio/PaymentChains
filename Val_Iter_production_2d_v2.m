%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Iteration 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize values
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
N_b=500; % grid points for debt level
N_k=400; % grid points for capital

% Value functions - Partial Equilibrium
beta= 0.8           ; % 0.8
R   = 1/beta        ;
q   = 1.25          ; % 1.75

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

% Special Values
b_bar=(f(K_low)-q*K_low)/(1-beta)    ; % one unit of labor input, or annuity value: natural borrowing limit
B_tilde=0.2*b_bar   ; % 0.4
Bstar=1/beta*(B_tilde-f(K_high)); % Threshold
Bast  = B_tilde/(q^(-1)-(q^(-1)-1)*(1-beta)*B_tilde);

% Bounds of Capital and Debt
scale=0.8;
K_lb =f_pinv(R*q*1.2);
K_ub =f_pinv(R*0.9);
B_min=Bstar-0.1;
B_max=b_bar*scale-1e-4;

% Vectorization
K_vec=linspace(K_lb,K_ub,N_k);
B_vec=linspace(B_min,B_max,N_b);

%% Functional Forms
% Consumption when all goods are bought at price 1
C_high = @(K,B) f(max(K(:),K_high)) - K_high - (B(:)' - B(:)' / R); % Upper bound is annuity plus put me at steady-state of capital
C_annuity_high=@(K,B) f(K(:)*0+K_high) - K_high - (B(:)' - B(:)' / R);

% Value function: no chaining
val_high = @(K,B) log(C_high(K,B))+beta*log(C_annuity_high(K,B)) / (1 - beta);

% Consumption when all goods are bought at price q
C_low = @(K,B) max((f(K(:)) - q * K_low - (B(:)' - B(:)' / R))/q,1e-6); % Upper bound is annuity plus put me at steady-state of capital
C_annuity_low=@(K,B) (f(K(:)*0+K_low) - q * K_low - (B(:)' - B(:)' / R))/q;

% Value function: all at price q
val_low = @(K,B) log(C_low(K,B))+beta*log(C_annuity_low(K,B)) / (1 - beta);

rango_plot=[K_lb K_ub B_min B_max];
if plotit==1
    figure('Name',"Value Function Bounds")
    fsurf(@(K,B) val_high(K,B),rango_plot,'FaceAlpha',0.5); hold on; axis tight;
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
% Pre-allocate storage
% valid = C_low(K_vec,B_vec) > 0 & isfinite(C_low(K,B));
V0 = val_low(K_vec,B_vec);
V_out = V0;
Bp_out = zeros(N_k, N_b);
Kp_out = zeros(N_k, N_b);
C_out = zeros(N_k, N_b);
S_out = zeros(N_k, N_b);
X_out = zeros(N_k, N_b);
val_gap = 2 * tol;

% Create all (Kp, Bp) combinations once
[Bp_mat, Kp_mat] = meshgrid(B_vec, K_vec); % size: [N_k x N_b]

% Parallel Processor
if isempty(gcp('nocreate'))
 %   parpool(4);  % or use parpool('local', 4);
end

%% Main Loop
while val_gap > tol
    V_new = zeros(N_k, N_b);
    C_out = V_new; S_out = V_new; X_out = V_new;
    Bp_out = V_new; Kp_out = V_new;

    parfor kk = 1:N_k
        K = K_vec(kk);
        for bb = 1:N_b
            B = B_vec(bb);

            % Vectorized over all Kp, Bp
            E_now = f(K) - B + Bp_mat / R;

            % Spot and chained goods
            S_now = min(max(B_tilde - B, 0), E_now);
            X_now = (E_now - S_now) / q;
            C_now = X_now + S_now - Kp_mat;

            % Mask valid consumption
            valid = imag(C_now) == 0 & C_now > 0 & isfinite(C_now);

            % Compute total value
            V_trial = -Inf(size(C_now));
            V_trial(valid) = log(C_now(valid)) + beta * V0(valid);

            % Max over valid entries
            [V_max, linIdx] = max(V_trial(:));
            V_new(kk, bb)   = V_max;

            % Extract indices
            [i_kp, i_bp] = ind2sub(size(V_trial), linIdx);

            % Store policies
            Bp_out(kk, bb) = B_vec(i_bp);
            Kp_out(kk, bb) = K_vec(i_kp);
            C_out(kk, bb)  = C_now(i_kp, i_bp);
            S_out(kk, bb)  = S_now(i_kp, i_bp);
            X_out(kk, bb)  = X_now(i_kp, i_bp);
        end
    end

    val_gap = max(abs(V_new(:) - V0(:))) / max(abs(V0(:)));
    V0 = V_new;
    if printinfo
        disp(['val gap: ', num2str(val_gap)]);
    end
end

%% Plotting Steady-State Solutions:
VH = val_high(K_vec, B_vec);  % size: [nK x nB]
VL = val_low(K_vec, B_vec);

% Mask out invalid (imaginary or Inf or NaN)
VH(~isreal(VH) | ~isfinite(VH)) = NaN;
VL(~isreal(VL) | ~isfinite(VL)) = NaN;

figure('Name','Value Function and Bounds');
surf(B_vec, K_vec, VH, 'FaceAlpha', 0.5, 'LineStyle', 'none'); hold on;
surf(B_vec, K_vec, VL, 'FaceAlpha', 0.5, 'LineStyle', 'none');
xlabel('Debt');
ylabel('Capital');
zlabel('Value');
title('Upper and Lower Bounds of Value Function');
axis tight; grid on;
surf(B_vec,K_vec,V_out,'FaceAlpha',0.5,'LineStyle',':'); drawnow; hold on;

figure('Name','Debt Prime')
surf(B_vec,K_vec',Bp_out,'FaceAlpha',0.5); drawnow; hold on;

figure('Name','Capital Prime')
surf(B_vec,K_vec',Kp_out,'FaceAlpha',0.5); drawnow; hold on;

figure;
quiver(B_vec,K_vec,Bp_out,Kp_out, 'k');

%% Plotting Steady-State Solutions:
figure('Name','Value Function')
surf(B_vec,K_vec',V_out,'FaceAlpha',0.5); drawnow; hold on;
xlabel('Debt'); ylabel('Capital'); 

figure('Name','Debt Prime')
surf(B_vec,K_vec',Bp_out,'FaceAlpha',0.5,'LineStyle','none'); drawnow; hold on;
xlabel('Debt'); ylabel('Capital'); 
surf(B_vec,K_vec',Bp_mat,'FaceAlpha',0.5,'LineStyle','none'); drawnow; hold on;

figure('Name','Debt Growth')
surf(B_vec,K_vec',Bp_out-Bp_mat,'FaceAlpha',0.5,'LineStyle','none'); drawnow; hold on;
xlabel('Debt'); ylabel('Capital'); 
surf(B_vec,K_vec',Bp_mat*0,'FaceAlpha',0.9,'LineStyle','none'); drawnow; hold on; 

figure('Name','Capital Prime')
surf(B_vec,K_vec',Kp_out,'FaceAlpha',0.5,'LineStyle','none'); drawnow; hold on;
xlabel('Debt'); ylabel('Capital'); 
surf(B_vec,K_vec',Kp_mat,'FaceAlpha',0.5,'LineStyle','none'); drawnow; hold on;

figure('Name','Capital Growth')
surf(B_vec,K_vec',Kp_out-Kp_mat,'FaceAlpha',0.5,'LineStyle','none'); drawnow; hold on;
xlabel('Debt'); ylabel('Capital'); 
surf(B_vec,K_vec',Kp_mat*0,'FaceAlpha',0.5,'LineStyle','none'); drawnow; hold on;

% Changes
dB = Bp_out - B_vec;     % Debt dynamics
dK = Kp_out - K_vec';    % Capital dynamics

figure;
contourf(B_vec,K_vec',dB,25); hold on;
step = 10;  % every 5th point
quiver(B_vec(1:step:end), K_vec(1:step:end), ...
       dB(1:step:end, 1:step:end), dK(1:step:end, 1:step:end),1.5, 'b');
xlabel('Debt'); ylabel('Capital'); 
title('Phase Diagram: Policy Dynamics');
axis tight; grid on;
line([B_tilde B_tilde],[K_vec(1) K_vec(end)],'LineWidth',2,'Color','r');
line([B_vec(1) B_vec(end)],[K_low K_low],'LineWidth',2,'Color','r','LineStyle','-');
line([B_vec(1) B_vec(end)],[K_high K_high],'LineWidth',2,'Color','r','LineStyle','-');
%line([Bstar Bstar],[K_vec(1) K_vec(end)],'LineWidth',2,'Color','r');

figure;
contourf(B_vec,K_vec',dK,25); hold on;
step = 10;  % every 5th point
quiver(B_vec(1:step:end), K_vec(1:step:end), ...
       dB(1:step:end, 1:step:end), dK(1:step:end, 1:step:end),1.5, 'b');
xlabel('Debt'); ylabel('Capital'); 
title('Phase Diagram: Policy Dynamics');
axis tight; grid on;
line([B_tilde B_tilde],[K_vec(1) K_vec(end)],'LineWidth',2,'Color','r');
line([B_vec(1) B_vec(end)],[K_low K_low],'LineWidth',2,'Color','r','LineStyle','-');
line([B_vec(1) B_vec(end)],[K_high K_high],'LineWidth',2,'Color','r','LineStyle','-');

%% Simulated Paths
% Simulation settings
T = 25;  % horizon

% Initial index positions (grid-aligned!)
init_ks = [N_k-10, 5, N_k - 10, 25;];
init_bs = [5, N_b-10, N_b - 10, N_b - 150];

% Plot setup
figure;
quiver(B_vec(1:step:end), K_vec(1:step:end), dB(1:step:end, 1:step:end), dK(1:step:end, 1:step:end), 0.5, 'k'); hold on;
xlabel('Debt'); ylabel('Capital');
title('Policy Vector Field with Trajectories');
grid on;

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
