%% main_paymentschain.m
% Main code for running key objects in payment-chain network model
clear; close all;

% Code Specs
printit=1; % Save plots into eps
plotit=1; % 1 if plot needed

%% Parameter Set
delta=0.9 ; % Value of discount in production
beta =0.95; % Discount factor
B_bar=1/(1-beta); 
fontsize=16;

% Colors
% Color Ordering for Plots
newcolors = [0.00 0.00 1.0
             0.5 0.5 0.9
             0.25 0.75 0.9
             0.25 0.20 0.9
             0.7 0.7 0.7];
color_t0=newcolors(1,:);
color_t1=newcolors(2,:);
color_t3=newcolors(3,:);

% Terminal Folder
params.folder='Matlab Figures\';

% Define base sizes
baseTick   = 18;      % tick labels
baseAxes   = 20;      % axes labels/title use multiplier
baseLegend = 20;      % legend text

% Set default fonts on all figures
set(groot, ...
    'DefaultLegendBox', 'off',...
    'DefaultAxesFontName','Times','DefaultAxesFontSize',baseTick, ...
    'DefaultAxesLabelFontSizeMultiplier',baseAxes/baseTick, ...
    'DefaultAxesTitleFontSizeMultiplier',baseAxes/baseTick, ...
    'DefaultLegendFontName','Times','DefaultLegendFontSize',baseLegend, ...
    'DefaultLegendFontSizeMode','manual',...
    'DefaultTextInterpreter','latex', ...
    'DefaultAxesTickLabelInterpreter','latex', ...
    'DefaultLegendInterpreter','latex', ...
    'DefaultLineLineWidth', 4);

%% Micro Objects
CalA_ex=@(mu,delta) delta*(1-mu)/(1-mu*delta);

%% Construction of Solutions
% mu vector
N_mu=1000;
min_mu=0;
max_mu=1;
mu_vec=linspace(min_mu,max_mu,N_mu);

% Two case comparison
N_mu_cases=3;
mu_cases=[0.7; 0.8; 0.85];
chain=0:1:8; % grid for distribution plots

% delta vector
N_delta=2;
delta_vec=[delta 0.5];

%% Vector of Solutions
Aeff_vec=zeros(N_delta,N_mu);
neff_vec=zeros(N_delta,N_mu);
peff_vec=zeros(N_delta,N_mu);
Apool_vec=zeros(N_delta,N_mu);
npool_vec=zeros(N_delta,N_mu);
ppool_vec=zeros(N_delta,N_mu);
A_vec=zeros(N_delta,N_mu);
for dd=1:N_delta
    for ii=1:length(mu_vec)
       delta_aux=delta_vec(dd);
       A_vec(dd,ii)=CalA_ex(mu_vec(ii),delta_aux);
       [Aeff_aux,p_aux,n_aux]=PlanCal2A_ex(mu_vec(ii),delta_aux);
       Aeff_vec(dd,ii)=Aeff_aux;
       neff_vec(dd,ii)=n_aux;
       peff_vec(dd,ii)=p_aux;
       [A_aux,p_aux,n_aux]=PoolCalA_ex(mu_vec(ii),delta_aux);
       Apool_vec(dd,ii)=A_aux;
       npool_vec(dd,ii)=n_aux;
       ppool_vec(dd,ii)=p_aux;
    end
end
rat_check=neff_vec(:,:)+(1-peff_vec(:,:));
eff_check=mu_vec./(1-mu_vec)-rat_check(1,:);

% Constructing Distribution - Efficient Case
mu1=mu_cases(1); mu2=mu_cases(2); mu3=mu_cases(3); 
pdf_2case=[(1-mu1)*mu1.^(chain); (1-mu2)*mu2.^(chain); (1-mu3)*mu3.^(chain)];

%% Efficient Assignment Plots
if plotit==1
    % Key micro objects - TFP Random Assignment (unternal use)
    figure('Name','TFP - Random Assignment')
    fplot(@(x) CalA_ex(x,0.9),[0.01 0.99],'LineWidth',2); hold on;
    fplot(@(x) CalA_ex(x,0.5),[0.01 0.99],'LineWidth',2);
    linestyleorder("mixedstyles")
    colororder(newcolors);
    xlabel('$\mu$','Interpreter','latex','FontSize',fontsize);
    ylabel('$\mathcal{A}(\mu)$','interpreter','latex','FontSize',fontsize); 
    grid on; axis tight; hold on;
    legend('High $\delta$','Low $\delta$','Interpreter','Latex','box','off')

    % Key micro objects - Positons Efficient Assignment (unternal use)
    figure('Name','Mininal Length - Random Assignment')
    plot(mu_vec,neff_vec,'LineWidth',2); hold on;
    linestyleorder("mixedstyles")
    colororder(newcolors);
    xlabel('$\mu$','Interpreter','latex');
    linestyleorder("mixedstyles")
    colororder(newcolors);

    % Main Comparison Figure - Probs Assignment (unternal use)
    figure('Name','Prob on minimal Length - Random Assignment')
    plot(mu_vec,peff_vec,'LineWidth',2); hold on; 
    linestyleorder("mixedstyles")
    colororder(newcolors);
    xlabel('$\mu$','Interpreter','latex');
    linestyleorder("mixedstyles")
    colororder(newcolors);

    % Main Comparison Figure
    figure('Name','TFP Comp')
    plot(mu_vec,A_vec,'LineWidth',2); hold on;
    plot(mu_vec,Aeff_vec,'LineWidth',2);
    xlabel('$\mu$','Interpreter','latex','FontSize',fontsize);
    ylabel('$\mathcal{A}(\mu,\delta)$','interpreter','latex','FontSize',fontsize); 
    grid on; axis tight; hold on;
    legend('Random $\mathcal{P}$ (high $\delta$)','Random $\mathcal{P}$  (low $\delta$)','Efficient $\mathcal{P}$  (high $\delta$)','Efficient $\mathcal{P}$  (low $\delta$)','Interpreter','Latex','box','off','location','southwest','AutoUpdate','Off')
    linestyleorder("mixedstyles")
    colororder(newcolors);
    ylim([0 1]); xlim([0 1]);
    line([0.5 0.5],[0 1],'LineWidth',0.5,'LineStyle','--'); hold on;
    if printit==1
        orient landscape;
        ax = gca;
        exportgraphics(ax,[params.folder 'TFP_example_comp.pdf'],'BackgroundColor','none');
    end

     % Main Comparison Figure
    figure('Name','TFP Comp + pool')
    plot(mu_vec,A_vec,'LineWidth',2); hold on;
    plot(mu_vec,Aeff_vec,'LineWidth',2);
    plot(mu_vec,Apool_vec,'LineWidth',2);
    xlabel('$\mu$','Interpreter','latex','FontSize',fontsize);
    ylabel('$\mathcal{A}(\mu)$','interpreter','latex','FontSize',fontsize); 
    grid on; axis tight; hold on;
    legend('random $\mathcal{P}$ (high $\delta$)','random $\mathcal{P}$  (low $\delta$)','efficient $\mathcal{P}$  (high $\delta$)','efficient $\mathcal{P}$  (low $\delta$)','pool $\mathcal{P}$  (high $\delta$)','pool $\mathcal{P}$  (low $\delta$)','Interpreter','Latex','box','off','location','southwest','AutoUpdate','Off')
    linestyleorder("mixedstyles")
    colororder(newcolors);
    ylim([0 1]); xlim([0 1]);
    line([0.5 0.5],[0 1],'LineWidth',0.5,'LineStyle','--'); hold on;

    % Key micro objects - Random Assignment
    figure('Name','TFP - Random Assignment')
    fplot(@(x) CalA_ex(x,0.9),[0.01 0.99],'LineWidth',2); hold on;
    fplot(@(x) CalA_ex(x,0.5),[0.01 0.99],'LineWidth',2); 
    xlabel('$\mu$','Interpreter','latex','FontSize',fontsize);
    ylabel('$\mathcal{A}(\mu)$','interpreter','latex','FontSize',fontsize); 
    grid on; axis tight; hold on;
    linestyleorder("mixedstyles")
    colororder(newcolors);
    legend('High $\delta$','Low $\delta$','Interpreter','Latex','box','off')
    if printit==1
        orient landscape;
        ax = gca;
     %   exportgraphics(ax,[params.folder 'TFP_example.pdf'],'BackgroundColor','none');
        % saveas(gcf, 'TFP_example', 'pdf')
    end

    figure('Name','Distributions - Random Assignment')       
    b=bar(chain,pdf_2case,'BarWidth',1.2); hold on;
    colororder(newcolors);
    xlabel('Chain Length','Interpreter','latex');
    ylabel('Probability Distribution','Interpreter','latex'); grid on; axis tight; 
    legend('low $\mu$','moderate $\mu$','high $\mu$','Interpreter','Latex','box','off');
     if printit==1
        orient landscape;
        % saveas(gcf, 'Dist_example', 'pdf')
        ax = gca;
        exportgraphics(ax,[params.folder 'DistRandom_example.pdf'],'BackgroundColor','none');

    end
end

% Distribution plot for efficient case
Aeff_2case=zeros(1,N_mu_cases);
neff_2case=zeros(1,N_mu_cases);
peff_2case=zeros(1,N_mu_cases);
pdfeff_2case=zeros(N_mu_cases,1)*chain;
Apool_2case=zeros(1,N_mu_cases);
npool_2case=zeros(1,N_mu_cases);
ppool_2case=zeros(1,N_mu_cases);
pdfpool_2case=zeros(N_mu_cases,1)*chain;
for mm=1:N_mu_cases
       delta_aux=delta_vec(dd);
       [Aeff_aux,p_aux,n_aux]=PlanCal2A_ex(mu_cases(mm),delta);
       Aeff_2case(mm)=Aeff_aux;
       neff_2case(mm)=n_aux;
       peff_2case(mm)=p_aux;
       pdfeff_2case(mm,n_aux+1)=p_aux; % adding one because support starts at zero in plots
       pdfeff_2case(mm,n_aux+2)=1-p_aux; % adding one because support starts at zero in plots
       [A_aux,p_aux,n_aux]=PoolCalA_ex(mu_cases(mm),delta);
       p_mix=1-(1-mu_cases(mm))^2 ;
       pdfpool_2case(mm,1)=1-p_mix;
       Apool_2case(mm)=A_aux;
       npool_2case(mm)=n_aux;
       ppool_2case(mm)=p_aux;
       pdfpool_2case(mm,n_aux+1)=p_aux*p_mix; % adding one because support starts at zero in plots
       pdfpool_2case(mm,n_aux+2)=(1-p_aux)*p_mix; % adding one because support starts at zero in plots
end


% Key micro objects - Efficient Assingment Assignment
if plotit==1
    figure('Name','TFP - Efficient Assignment')
    plot(mu_vec,Aeff_vec,'LineWidth',2); hold on;
    xlabel('$\mu$','Interpreter','latex','FontSize',fontsize);
    ylabel('$\mathcal{A}(\mu)$','interpreter','latex','FontSize',fontsize); 
    grid on; axis tight; hold on;
    legend('High $\delta$','Low $\delta$','Interpreter','Latex','box','off')
    if printit==1
        orient landscape;
        ax = gca;
        exportgraphics(ax,[params.folder 'TFP_example.pdf'],'BackgroundColor','none');       
    end

    figure('Name','Distributions - Random Assignment')
    b=bar(chain,pdfeff_2case,'BarWidth',1.2); hold on;       
    colororder(newcolors);
    xlabel('Chain Length','Interpreter','latex');
    ylabel('Probability Distribution','Interpreter','latex'); grid on; axis tight; 
    legend('low $\mu$','moderate $\mu$','high $\mu$','Interpreter','Latex','box','off');
     if printit==1
        orient landscape;
        ax = gca;
        exportgraphics(ax,[params.folder 'DistEff_example.pdf'],'BackgroundColor','none');
     end

    % Pooled case
    figure('Name','Distributions - Pool Assignment')
    b=bar(chain,pdfpool_2case,'BarWidth',1.2); hold on;       
    colororder(newcolors);
    xlabel('Chain Length','Interpreter','latex');
    ylabel('Probability Distribution','Interpreter','latex'); grid on; axis tight; 
    legend('low $\mu$','moderate $\mu$','high $\mu$','Interpreter','Latex','box','off');
     if printit==1
        orient landscape;
        ax = gca;
        exportgraphics(ax,[params.folder 'DistPool_example.pdf'],'BackgroundColor','none');
    end
end


%% Functions 
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

function [A,p,n]=PlanCal2A_ex(mu,delta)
    n=floor(mu/(1-mu));
    n1=ceil(mu/(1-mu));
    if n==0
        A=delta;
        p=1-mu/(1-mu);
        n=0;
    else        
        p=n1-mu/(1-mu);
        w_short=p/(p*n+(1-p)*n1);
        w_long=(1-p)/(p*n+(1-p)*n1);
        A_short=delta*(1-delta^(n))/(1-delta)*w_short;
        A_long=delta*(1-delta^(n1))/(1-delta)*w_long;
        A=A_short+A_long;
    end
end

function [A,p,n]=PoolCalA_ex(mu,delta)
    rat=1/(1-mu);
    n=floor(rat);
    n1=ceil(rat);
    if n==0
        A=delta;
        p=1-rat;
        n=0;
    else        
        p=n1-rat;
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