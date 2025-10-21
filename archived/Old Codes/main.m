%% main_paymentschain.m
% Main code for running key objects in payment-chain network model

%% Micro Objects
params.folder='Matlab Figures\';
CalA_ex=@(mu,delta) delta*(1-mu)/(1-delta*mu);

% Key micro objects
if plotit==1
    figure
    fplot(@(x) CalA_ex(x,0.9),[0.01 0.99],'LineWidth',2,'Color',color_t0); hold on;
    fplot(@(x) CalA_ex(x,0.5),[0.01 0.99],'LineWidth',2,'Color',color_t1); xlabel('$\mu$','Interpreter','latex','FontSize',fontsize);
    ylabel('$\mathcal{A}(\mu)$','interpreter','latex','FontSize',fontsize); 
    grid on; axis tight; hold on;
    legend('High $\delta$','Low $\delta$','Interpreter','Latex','box','off')
    if printit==1
        orient landscape;
        ax = gca;
        exportgraphics(ax,[params.folder 'TFP_example.pdf'],'BackgroundColor','none');
        % saveas(gcf, 'TFP_example', 'pdf')
    end

    figure;
    chain=0:1:8;
    mu1=0.6; mu2=0.7;
    b=bar(chain,[(1-mu1)*mu1.^(chain);(1-mu2)*mu2.^(chain)],'BarWidth',1.2); hold on;
    b(1).FaceColor=color_t0;
    b(2).FaceColor=color_t1; xlabel('Chain Length','Interpreter','latex');
    ylabel('Probability Distribution','Interpreter','latex'); grid on; axis tight; 
    legend('Low $\mu$','High $\mu$','Interpreter','Latex','box','off');
     if printit==1
        orient landscape;
        % saveas(gcf, 'Dist_example', 'pdf')
        ax = gca;
        exportgraphics(ax,[params.folder 'Dist_example.pdf'],'BackgroundColor','none');

    end
end

q_mu=@(mu) CalA(mu,0.9)^(-1);
if plotit==1
    figure
    fplot(@(mu) q_mu(mu), [0.01 0.99]); xlabel('$\mu$','Interpreter','latex');
    title('$q(\mu)$','Interpreter','latex'); grid on; axis tight;
end

% mu vector
mu_vec=linspace(0.01,0.9999,1000);
delta=0.9;
for ii=1:length(mu_vec);
   A_vec(ii)=CalA_man(mu_vec(ii),delta);
   Aeff_vec(ii)=PlanCal2A_ex(mu_vec(ii),delta);
   Aplan_vec(ii)=ExtPlanCalA_ex(mu_vec(ii),delta);
end

if plotit==1
    figure('Name','Comp')
    plot(mu_vec,A_vec,'LineWidth',2,'Color',color_t0); hold on;
    plot(mu_vec,Aplan_vec,'LineWidth',2,'Color',color_t1);
    plot(mu_vec,Aeff_vec,'LineWidth',1,'Color','r');
    xlabel('$\mu$','Interpreter','latex','FontSize',fontsize);
    ylabel('$\mathcal{A}(\mu)$','interpreter','latex','FontSize',fontsize); 
    grid on; axis tight; hold on;
    legend('Baseline','Efficient','Interpreter','Latex','box','off')

end

%% Plotting Efficient Solution
function [A]=EffCalA_ex(mu,delta)
    n=floor(mu/(1-mu));
    r=mod(mu/(1-mu),1);
    if n==0
        A=delta;
    else        
        cum=0;
        for i=1:n
            cum=cum+(delta)^(i);
        end
        A=(n*(1/n*cum)+(delta)^(n+1)*r)/(mu/(1-mu));
    end
end

%% Plotting Efficient Solution
function [A]=PlanCalA_ex(mu,delta)
    n=floor(mu/(1-mu));
    n1=ceil(mu/(1-mu));
    if n==0
        A=delta;
    else        
        p=n1-mu/(1-mu);
        cum=0;
        for i=1:n
            cum=cum+(delta)^(i);
        end
        A_short=cum/n;
        A_long=(cum+delta^(n1))/n1;
        A=A_short*p+A_long*(1-p);
    end
end

function [A]=PlanCal2A_ex(mu,delta)
    n=floor(mu/(1-mu));
    n1=ceil(mu/(1-mu));
    if n==0
        A=delta;
    else        
        p=n1-mu/(1-mu);
        w_short=p/(p*n+(1-p)*n1);
        w_long=(1-p)/(p*n+(1-p)*n1);
        A_short=delta*(1-delta^(n))/(1-delta)*w_short;
        A_long=delta*(1-delta^(n1))/(1-delta)*w_long;
        A=A_short+A_long;
    end
end

function [A]=PlanCal3A_ex(mu,delta)
        A=delta*(1-delta^(max(mu/(1-mu),1)))/(1-delta)/max(mu/(1-mu),1);
end

function [A]=ExtPlanCalA_ex(mu,delta)
    G1=min((1-mu)/mu,1);
    Gn=max(1-G1,0);
    if Gn==0 
        n =0;
        A=delta;
    else        
        n=((mu)/(1-mu)-G1)/Gn;
        A=delta*G1+delta*(1-delta^(n))/(1-delta)*Gn/n;
    end
end

function [A]=CalA_man(mu,delta)
        A=(1-mu)*delta/(1-mu*delta);
end