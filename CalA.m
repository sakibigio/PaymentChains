function A_out=CalA(mu,delta) 
    if mu==1
        A_out=10e-16;
    elseif mu==0
        A_out=delta;
    else
        A_out=delta/(1-delta)*(1-mu)./mu.*log((1-mu.*delta)./(1-mu));
    end
end