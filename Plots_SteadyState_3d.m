%% Steady-State Plots
B_tilde= params.B_bar/5;
B_vec  = linspace(Bstar(B_tilde)-1,Bast(B_tilde)+2,2000)  ;
R_vec  = NaN*B_vec; B_p_vec  = NaN*B_vec; val_vec  = NaN*B_vec; flag_vec  = NaN*B_vec;
ii=0;
for B0=B_vec
    ii=ii+1;
    [B_p_vec(ii),val_vec(ii),flag_vec(ii)]=B_prima(B0,B_tilde,B_tilde,funcs);
    R_vec(ii)=1/params.beta*B_p_vec(ii)/B0;
end

figure('Name','B 45-degree plot')
plot(B_vec,B_p_vec,'LineWidth',2); grid on; axis tight; hold on;
plot(B_vec,B_vec,'LineWidth',2,'LineStyle','--'); grid on; axis tight;
grid on; yvals=ylim; line([Bstar(B_tilde) Bstar(B_tilde)],[yvals(1) yvals(2)],'Color','k','LineWidth',2,'LineStyle','--');
grid on; yvals=ylim; line([Bast(B_tilde_p) Bast(B_tilde_p)],[yvals(1) yvals(2)],'Color','k','LineWidth',2,'LineStyle','--');
grid on; yvals=ylim; line([B_tilde B_tilde],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.5],'LineWidth',2,'LineStyle','-.');
grid on; yvals=ylim; line([Bhat(B_tilde,B_tilde) Bhat(B_tilde,B_tilde)],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.5],'LineWidth',2,'LineStyle','-.');
grid on; yvals=ylim; line([Bhat(B_tilde,B_tilde) Bhat(B_tilde,B_tilde)],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.5],'LineWidth',2,'LineStyle','-.');
xlabel('B'); ylabel('B''');

figure('Name','Rate Plots')
index1=(B_vec<Bstar(B_tilde));
index2=((B_vec>Bstar(B_tilde))&(B_vec<Bast(B_tilde)));
index3=(B_vec>Bast(B_tilde));
B_vec<Bast(B_tilde)
plot(B_vec(index1),R_vec(index1),'LineWidth',2,'Color',color_t0); axis tight; grid on; hold on;
plot(B_vec(index2),R_vec(index2),'LineWidth',2,'Color',color_t0); grid on;  
plot(B_vec(index3),R_vec(index3),'LineWidth',2,'Color',color_t0); grid on;  
% title('$R(B)$','interpreter', 'latex'); 
% legend('$\mathcal{E}_{t}$','$\mathcal{E}_{t+1}$','Box','Off','interpreter', 'latex','Location','northwest','AutoUpdate','off')
grid on; yvals=ylim; line([Bstar(B_tilde) Bstar(B_tilde)],[yvals(1) yvals(2)],'Color','k','LineWidth',2,'LineStyle','--');
grid on; yvals=ylim; line([Bast(B_tilde) Bast(B_tilde)],[yvals(1) yvals(2)],'Color','k','LineWidth',2,'LineStyle','--');
grid on; yvals=ylim; line([B_tilde B_tilde],[yvals(1) yvals(2)],'Color',[0.5 0.5 0.5],'LineWidth',2,'LineStyle','-.');
grid on; yvals=ylim; xvals=xlim; line([xvals(1) xvals(end)],[1/params.beta 1/params.beta],'Color',[0.3 0.3 0.3],'LineWidth',2,'LineStyle',':');
grid on; yvals=ylim; line([Bhat(B_tilde,B_tilde) Bhat(B_tilde,B_tilde)],[yvals(1) yvals(2)],'Color',[0.3 0.3 0.3],'LineWidth',2,'LineStyle','-.');
B_star_aux=Bstar(B_tilde);
set(gca,'XTick',B_star_aux,'XTickLabel',[]);
set(gca,'YTick',[yvals(1) yvals(2)],'YTickLabel',[]);
yTicks = get(gca,'ytick');
xTicks = get(gca, 'xtick');

ax = axis; %Get left most x-position
text(ax(1)-0.01,yTicks(2),'$\beta^{-1}$', 'HorizontalAlignment','Right','interpreter', 'latex');
text(xTicks(1),yTicks(1)-0.005,'$B^{\star}$', 'HorizontalAlignment','Center','VerticalAlignment','Middle','interpreter', 'latex');
text(Bast(B_tilde),yTicks(1)-0.005,'$B^{\ast}(\tilde{B})$', 'HorizontalAlignment','Center','VerticalAlignment','Middle','interpreter', 'latex');
text(B_tilde,yTicks(1)-0.005,'$\tilde{B}$', 'HorizontalAlignment','Center','VerticalAlignment','Middle','interpreter', 'latex'); hold on;
if printit==1
    orient landscape;
   %  saveas(gcf,'F_SSRate','pdf');
   ax = gca;
    exportgraphics(ax,[params.folder 'F_SSRate.pdf'],'BackgroundColor','none');
end

%% Plot for R
% fix arbitrary debt level. 
B0=params.B_bar/5;
B_vec  = linspace(0.5,1+params.beta*B0+1,100)  ;
Bp_vec = linspace(0.5,1+params.beta*B0+1.1,50)   ;
B_p_mat = NaN(length(B_vec),length(Bp_vec));
R_mat   = NaN(length(B_vec),length(Bp_vec));
val_mat=zeros(length(B_vec),length(Bp_vec))  ;
flag_mat=zeros(length(B_vec),length(Bp_vec)) ;
%gamma_mat=zeros(length(B_vec),length(Bp_vec));
%rho_mat=zeros(length(B_vec),length(Bp_vec))  ;
%eta_mat=zeros(length(B_vec),length(Bp_vec))  ;
%one_mat=zeros(length(B_vec),length(Bp_vec))  ;


% Some Tests
%mu_1=mu_Z(B0,1+beta*B0);
%q_1 =q_Z(B0,1+beta*B0);
%mu_test = mu_Z(B0,B_vec); 
%ii=0; q_test=zeros(length(mu_test),1); q_test2=zeros(length(mu_test),1);
%for B_tilde=B_vec
%    ii=1+ii;
%    q_test(ii)  = q_Z(B0,B_tilde);
%    q_test2(ii) = (q_Z(B0,B_tilde)-1)*(B_tilde-B0);
%end

% Main Loop
ii=0;
%R_iter=1/beta;
for B_tilde=B_vec
    jj=0;
    ii=ii+1;
    for B_tilde_p=Bp_vec
        jj=jj+1;
        [B_p_mat(ii,jj),val_mat(ii,jj),flag_mat(ii,jj)]=B_prima(B0,B_tilde,B_tilde_p,funcs);
        R_mat(ii,jj)=1/params.beta*B_p_mat(ii,jj)/B0;
    end
end

%% Key Figures - Phase Diagram
figure('Name','3d Rate plot')
surf(B_vec,Bp_vec,params.beta*R_mat');
title('$\beta \cdot R(B,\tilde{B},\tilde{B}^{''})$','Interpreter','latex');
xlabel('$\tilde{B}$','Interpreter','latex');
ylabel('$\tilde{B}^{''}$','Interpreter','latex'); axis tight;

figure('Name','3d Rate plot')
imagesc(B_vec,Bp_vec,R_mat');
title('$\beta \cdot R(B,\tilde{B},\tilde{B}^{''})$','Interpreter','latex');
xlabel('$\tilde{B}$','Interpreter','latex');
ylabel('$\tilde{B}^{''}$','Interpreter','latex'); axis tight;

figure('Name','error plot')
surf(B_vec,Bp_vec,val_mat');
title('Numerical Error','Interpreter','latex');
xlabel('$\tilde{B}$','Interpreter','latex');
ylabel('$\tilde{B}^{''}$','Interpreter','latex');


figure
surf(B_vec,Bp_vec,flag_mat');
title('Flags','Interpreter','latex');
xlabel('$\tilde{B}$','Interpreter','latex');
ylabel('$\tilde{B}^{''}$','Interpreter','latex');