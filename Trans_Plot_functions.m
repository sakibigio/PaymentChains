%% Equilibrium Rate Function
function [R_out,fval,flag]=R(B,B_tilde,B_tilde_p,params)
    beta=params.beta;
    e_cond=(B_tilde>B&&B_tilde<1+beta*B); % expenditure condition
    s_cond=(B_tilde_p>B&&B_tilde_p<1+beta*B); % expenditure condition
    if (e_cond==1)&&(s_cond==1)
        [x,fval,flag]=gamma(B,B_tilde,B_tilde_p,params);
        R_out=1/beta*x;
    elseif (e_cond==1)&&(s_cond==0)
        [x,fval,flag]=eta(B,B_tilde,params);
        R_out=1/beta*x;
    elseif (e_cond==0)&&(s_cond==1)
        [x,fval,flag]=rho(B,B_tilde_p,params);
        R_out=1/beta*x;
    else 
        R_out=1/beta;
        fval=0; flag=1;
    end
end

function [Rg_out,fval,flag]=Rg(B,B_tilde,B_tilde_p,params)
    beta=params.beta;
    [Gamma,fval,flag]=gammaG(B,B_tilde,B_tilde_p,params);
    Rg_out=1/beta*Gamma;
end



