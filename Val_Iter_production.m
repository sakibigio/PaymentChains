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
% Value functions - Partial Equilibrium
beta= 0.8           ; % 0.8
R   = 1/beta        ;
q   = 1.75          ; % 1.75

% Production Function
alpha=0.5; % DRS
f=@(k) k^(alpha)  ;
f_inv=@(x) x^(1/alpha);
f_p=@(k) alpha*k^(alpha-1); % f prime
f_pinv=@(f) (f/alpha)^(1/(alpha-1)); % f prime

% We discussed two possible models, one where agents buy capital in
% the market and one where they don't
% Version: capital bought in the market:
K_high=f_pinv(R)  ;
K_low =f_pinv(R*q);

% Special Forms
b_bar=(f(K_low)-q*K_low)*beta/(1-beta)    ; % one unit of labor input, or annuity value: natural borrowing limit
B_tilde=0.4*b_bar   ; % 0.4
Bstar=1/beta*(B_tilde-1); % Threshold
Bast  = B_tilde/(q^(-1)-(q^(-1)-1)*(1-beta)*B_tilde);

% Support of values of debt
a_max=f(K_high)-q*K_high;
a_bar=f(K_low)-R*b_bar;
range=[a_bar+(10e-10) a_max*2]    ;
A_vec=linspace(range(1),range(end),1000); % grid space

%% Functional Forms
% Extreme values
C_high=@(A)   (A-K_high)+b_bar; % annuity value of consumption
val_high=@(A) log(C_high(A))/(1-beta)  ; % no chained consumption
C_low=@(A)    (A+b_bar-q*K_low)/q; % annuity value of consumption
val_low=@(A)  log(C_low(A))/(1-beta) ; % all consumption is at chained price q

if plotit==1
    figure('Name',"Value Function Bounds")
    fplot(@(A) val_high(A),range,'LineWidth',3,'LineStyle','-'); hold on; axis tight;
    fplot(@(A) val_low(A),range,'LineWidth',3,'LineStyle','--'); grid on;% 
    figure('Name',"Consumption")
    fplot(@(A) C_high(A),range,'LineWidth',3,'LineStyle','-'); hold on; axis tight;
    fplot(@(A) C_low(A),range,'LineWidth',3,'LineStyle','--'); grid on;% 
end


%% Value Function Iteration
% define functions - B_p stands for B_prime (future debt), determining expenditures today
E     = @(B_p,A) A+B_p                    ; % total expenditures as functions of current wealth minus debt

% knowing expenditures and credit line, you can determine optimal
% expenditures in spot and chained, and determine total consumption
S_w   = @(E,A,B_tilde) min([max(B_tilde+A,0); E]);
X_w   = @(E,A,B_tilde) (E-S_w(E,A,B_tilde))/q ;
C_w   = @(E,A,B_tilde,i) X_w(E,A,B_tilde)+S_w(E,A,B_tilde)-i; % Consumption is goods purchased minus i investment

% Main Iteration 
% initialize values
V0=val_low(A_vec); V_out=V0; C_out=0*V_out; S_out=0*V_out; X_out=0*V_out; Bp_out=0*V_out; i_out=0*V_out;
v_out_ingrid=0*V_out; bp_ingrid=0*V_out;
% %% A_prime on grid 
val_gap=2*tol; % reset tolerance
while val_gap>tol
    if printinfo
        plot(A_vec,V0,'LineWidth',1,'LineStyle',':'); drawnow;
        disp(['Val at iter:' num2str(val_gap)]);
    end
    aa=0;
    for A=A_vec
        aa=aa+1;

        % Update definition of RHS of Value Functions
        pp=0;
        for a_p=A_vec
            pp=pp+1;
            feasible=(a_p+R*b_bar>0)&&(C_w(E(b_bar,A),A,B_tilde,f_inv(a_p+R*b_bar))>0);
            if feasible
                Tv=@(bp) -1*(log(C_w(E(bp,A),A,B_tilde,f_inv(a_p+R*bp)))+beta*V0(pp));
    
                % Solution
                bp_g=b_bar;
                % [X_aux,v_out_aux]=fminsearch(Tv,bp_g);
                [X_aux,v_out_aux]=fmincon(Tv,bp_g,[],[],[],[],-a_p/R,b_bar);
                v_out_ingrid(pp)=-v_out_aux;
                bp_ingrid(pp)=X_aux;
            else
                v_out_ingrid(pp)=1i;
            end
        end
        Tv_iter=v_out_ingrid;
        index=(imag(Tv_iter)==0);

        % evaluate at optimum to obtain solution
        [V_out(aa),I]=max(Tv_iter(index));
        Bp_aux=bp_ingrid(I);

        % update solutions
        Bp_out(aa)= Bp_aux;
        i_out(aa) = f_inv(a_p+R*Bp_aux);
        C_out(aa) = C_w(E(Bp_out(aa),A),A,B_tilde,i_out(aa));
        S_out(aa) = S_w(E(Bp_out(aa),A),A,B_tilde);
        X_out(aa) = X_w(E(Bp_out(aa),A),A,B_tilde);
    end
    val_gap=abs(max(V_out-V0));
    V0=V_out;
end
% Numerical Exit point
B_exit_num=min(B_vec(abs(V0-val_low(B_vec))<10e-10));


%% Double Optimization
val_gap=2*tol; % reset tolerance
while val_gap>tol
    if printinfo
        plot(A_vec,V0,'LineWidth',1,'LineStyle',':'); drawnow;
        disp(['Val at iter:' num2str(val_gap)]);
    end
    aa=0;
    for A=A_vec
        aa=aa+1;
        % Update definition of RHS of Value Functions
        Tv=@(X) -1*(log(C_w(E(X(1),A),A,B_tilde,exp(X(2))))+beta*interp1(A_vec,V0,max(f(exp(X(2)))-X(1),range(1)),'nearest'));
        
        % Solution
        guess=[-(A_vec(1)-f(K_low)),log(K_low)];
        [X_aux,v_out_aux]=fminsearch(Tv,guess);

        % update solutions
        Bp_out(aa)= X_aux(1);
        i_out(aa) = exp(X_aux(2));
        V_out(aa) = -v_out_aux;
        C_out(aa) = C_w(E(Bp_out(aa),A),A,B_tilde,i_out(aa));
        S_out(aa) = S_w(E(Bp_out(aa),A),A,B_tilde);
        X_out(aa) = X_w(E(Bp_out(aa),A),A,B_tilde);
    end
    val_gap=abs(max(V_out-V0));
    V0=V_out;
end
% Numerical Exit point
B_exit_num=min(B_vec(abs(V0-val_low(B_vec))<10e-10));

