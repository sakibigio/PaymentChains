function q_out=q_mu(mu,delta) 
    T=length(mu); q_out=zeros(T,1);
    for tt=1:T
        if mu(tt)==1
            q_out(tt)=10e16;
        elseif mu(tt)==0
            q_out(tt)=1/delta;
        elseif mu(tt)<10e-10
            q_out(tt)=delta;
        else
            q_out(tt)=CalA(mu(tt),delta).^(-1);
        end
    end
end