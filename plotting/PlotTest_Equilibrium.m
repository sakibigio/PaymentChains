%% Key Plots for Equilibrium objects
B0=1.5;
range_test=[0.9 3]; 

B_tilde=2; B_tilde_p=2;
figure('Name','Test Marginal Price')
fplot(@(x) q_E(x,B_tilde),range_test,'LineWidth',2); title('Borrowing Marginal Price'); axis tight;
grid on; yvals=ylim; line([B_tilde B_tilde],[yvals(1) yvals(2)],'Color','k','LineWidth',2);

figure
fplot(@(x) q_B(x,B_tilde),range_test,'LineWidth',2); title('Borrowing Marginal Price');
grid on; yvals=ylim; line([B_tilde B_tilde],[yvals(1) yvals(2)],'Color','k','LineWidth',2);

figure
fplot(@(x) Q_Z(x,B_tilde_p),range_test,'LineWidth',2); title('Average Price');
grid on; yvals=ylim; line([B_tilde B_tilde],[yvals(1) yvals(2)],'Color','k','LineWidth',2);

figure
fplot(@(x) Q_Z(x,B_tilde_p)./Q_Z(B0,B_tilde),range_test,'LineWidth',2); title('Ratio of Average Prices Ratio');
grid on; yvals=ylim; line([B_tilde B_tilde],[yvals(1) yvals(2)],'Color','k','LineWidth',2);

figure
fplot(@(x) Pi(2,x,B_tilde,B_tilde_p),range_test,'LineWidth',2); title('Marginal Inflation');
grid on; yvals=ylim; line([B_tilde B_tilde],[yvals(1) yvals(2)],'Color','k','LineWidth',2);
