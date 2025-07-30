function [hist]=Val_iter_sim(T,init_bs,B_vec,policies)
% Unpack policies
C_out=policies.C_out;
S_out=policies.S_out;
X_out=policies.X_out;
Bp_out=policies.Bp_out;

% Number of cases
N_c=length(init_bs);

%  Set time
time=1:1:T;

% Clear History
B_hist=zeros(T,N_c); C_hist=B_hist; S_hist=B_hist; X_hist=B_hist;

% Main simulation
    for i = 1:N_c
        b_idx = init_bs(i);
    
        B_path = zeros(T, 1);
        C_path = zeros(T, 1);
        S_path = zeros(T, 1);
        X_path = zeros(T, 1);
        for t = time
            B_path(t) = B_vec(b_idx);
           
            % Advance to next state using grid indices
            B_next = Bp_out(b_idx);
    
            % Update Consumption
            C_path(t) = C_out(b_idx);
            X_path(t) = X_out(b_idx);
            S_path(t) = S_out(b_idx);
    
            % Find indices for next step
            [~, b_idx] = min(abs(B_vec - B_next));
        end
    
       % Update histories
       B_hist(:,i)=B_path;
       C_hist(:,i)=C_path;
       X_hist(:,i)=X_path;
       S_hist(:,i)=S_path;
    end
hist.B_hist=B_hist;
hist.C_hist=C_hist;
hist.X_hist=X_hist;
hist.S_hist=S_hist;
end