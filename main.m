%% main_paymentschain.m
% Main code for running key objects in payment-chain network model

%% Micro Objects
% Key micro objects
if plotit==1
    CalA_ex=@(mu,delta) delta/(1-delta)*(1-mu)./mu.*log((1-mu.*delta)./(1-mu));
    figure
    fplot(@(x) CalA_ex(x,0.9),[0.01 0.99],'LineWidth',2,'Color',color_t0); hold on;
    fplot(@(x) CalA_ex(x,0.5),[0.01 0.99],'LineWidth',2,'Color',color_t1); xlabel('$\mu$','Interpreter','latex','FontSize',fontsize);
    title('$\mathcal{A}(\mu)$','interpreter','latex','FontSize',fontsize); grid on; axis tight; hold on;
    legend('High $\delta$','Low $\delta$','Interpreter','Latex')
    if printit==1
        orient landscape;
        ax = gca;
        exportgraphics(ax,'TFP_example.pdf','BackgroundColor','none');
        % saveas(gcf, 'TFP_example', 'pdf')
    end

    figure;
    chain=0:1:8;
    mu1=0.6; mu2=0.7;
    b=bar(chain,[(1-mu1)*mu1.^(chain);(1-mu2)*mu2.^(chain)],'BarWidth',1.2); hold on;
    b(1).FaceColor=color_t0;
    b(2).FaceColor=color_t1; xlabel('Chain Length','Interpreter','latex');
    title('Probability Distribution','Interpreter','latex'); grid on; axis tight; 
    legend('Low $\mu$','High $\mu$','Interpreter','Latex');
     if printit==1
        orient landscape;
        % saveas(gcf, 'Dist_example', 'pdf')
        ax = gca;
        exportgraphics(ax,'Dist_example.pdf','BackgroundColor','none');

    end
end

% q_mu=@(mu) CalA(mu);
if plotit==1
    figure
    fplot(q_mu,[0.01 0.99]); xlabel('$\mu$','Interpreter','latex');
    title('$q(\mu)$','Interpreter','latex'); grid on; axis tight;
end

%% Computing and Plotting Value Functions

%% Computing and Plotting Transitions and Comparison with Planner Problem

%% Analysis of Resilience