% CT_analyticplots
figure('Name','b_h plot')
fplot(@(x) RHS_b_h(x),[-30,30]); hold on
fplot(@(x) LHS_b_h(x),[-30,30]); hold on

% Plotting Value function
figure('Name','Value Function (Continuous-Time Solution)')
hold on;
plot(s_vec(index_aux),URF((w1+rho*b_h)/q)/rho+v_bar*(s_vec(index_aux)-b_h),'LineWidth',2,'LineStyle','-','Color','k');
plot(s_vec(index_aux2),URF((w1+rho*s_vec(index_aux2))/q)/rho,'LineWidth',2,'LineStyle','-','Color','k');
plot(s_vec(s_unc_index),URF((w1+rho*s_vec(s_unc_index)))/rho,'LineWidth',2,'LineStyle','-','Color','k');
grid on; axis tight;

figure('Name','Consumption Function (Continuous-Time Solution)')
plot(s_vec(index_aux),c_back+0*s_vec(index_aux),'Color','b'); hold on; grid on;
plot(s_vec(index_aux2),(w1+rho*s_vec(index_aux2))/q,'Color','b'); hold on; grid on;
plot(s_vec(s_unc_index),w1+rho*s_vec(s_unc_index),'Color','b'); 
line([s_vec(1) s_vec(end)],[c_back c_back],'LineStyle',':','LineWidth',1); 
scatter(s_bl,c_star,40,'MarkerFaceColor','b')
scatter(s_bl,c_back,40,'MarkerFaceColor','w')
scatter(b_h,(w1+rho*b_h)/q,40,'MarkerFaceColor','b')
scatter(b_h,c_back,40,'MarkerFaceColor','w')
ylimits=ylim;
line([b_h b_h],[ylimits(1) ylimits(2)],'LineWidth',1);
line([s_bl s_bl],[ylimits(1) ylimits(2)],'LineWidth',1);
axis tight;