% CT_solution
c_star=w1+rs*s_bl;
c_back=fsolve(@(c) Uprime(c)/q-(URF(c_star)-URF(c))/(w1-q*c+rho*s_bl),w1/25);
c_drop=c_back;
v_bar = (URF(c_star)-URF(c_back))/(w1-q*c_back+rho*s_bl);

% Solve b_h
RHS_b_h=@(b_h) URF((w1+rho*b_h)/q);
LHS_b_h=@(b_h) URF(c_back)+v_bar*(w1-q*c_back+rho*b_h);
b_h=fsolve(@(x) RHS_b_h(x)-LHS_b_h(x),s_bl);
index_aux=((s_vec>b_h)&(s_vec<s_bl));
index_aux2=((s_vec<=b_h));