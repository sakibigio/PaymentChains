function sim=sim_planner(B_tilde_t,theta,params,funcs,B0)
    % Parameters
    beta=params.beta; 
    delta=params.delta;

    % Begin main loop
    y_check=funcs.y_check;
    income_check=funcs.income_check;
    C_s=funcs.C_s;
    S_w=funcs.S_w;
    mu_Z=funcs.mu_Z;
 %   q_Z=funcs.q_Z;
 %   X_w=funcs.X_w;
    Y_Z=funcs.Y_Z;
    P=funcs.P;
    B_Ram=funcs.B_Ram;

    % Begin main simulation
    T=length(B_tilde_t);
    B_t=B0*ones(T,1);
    Rg_t=ones(T,1); valg_t=zeros(T,1); flagg_t=zeros(T,1); 
    for tt=1:T
        if tt==1
            B0=B_t(tt); B_tilde=B_tilde_t(tt); B_tilde_p=B_tilde_t(tt+1); % Update State Variable
            [B_t(tt+1),valg_t(tt),flagg_t(tt)]=B_prima(B0,B_tilde,B_tilde_p,funcs);
            Rg_t(tt)=B_t(tt+1)/B_t(tt); 
        else
             B_tilde=B_tilde_t(tt); 
            B_t(tt)=B_Ram(B_tilde,theta);
        end
    end   
    C_s_t = C_s(B_t);
    S_w_t = S_w(B_t,B_tilde_t);
    mu_t  = mu_Z(B_t,B_tilde_t);
    q_t   = q_mu(mu_t,delta);
    X_w_t  = mu_t./q_t;
    Y_t   = Y_Z(B_t,B_tilde_t);
    y_check_t=y_check(B_t,B_tilde_t) ;
    income_check_t=income_check(B_t,B_tilde_t);

    % Probably we want to indicate the regions too

    % Saver Variables
    sim.B_t=B_t;
    sim.B_tilde_t=B_tilde_t;
    sim.R_t=1/beta*Rg_t;
    sim.C_s_t=C_s_t;
    sim.S_w_t=S_w_t;
    sim.mu_t=mu_t;
    sim.q_t=q_t;
    sim.X_w_t=X_w_t;
    sim.Y_t=Y_t;    
    sim.y_check_t=y_check_t;
    sim.income_check_t=income_check_t;
end