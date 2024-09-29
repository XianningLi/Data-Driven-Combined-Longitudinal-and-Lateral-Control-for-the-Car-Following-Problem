State_tra = load('State_tra.csv');
load('Su.mat','Su');
T = 0:0.01:size(State_tra,1)*0.01-0.01;
State_tra_com = load('State_tra_com.csv');

%% error of the ADP controller
Q = 10*eye(6);
R = 1*eye(3);
Cost = 0;
err = zeros(length(State_tra),6);
for i=1:length(State_tra)
    err(i,:) = get_error(State_tra(i,:))';
    ue = Su*get_basis_ue(State_tra(i,:));
    Cost = Cost + (ue'*R*ue + err(i,:)*Q*err(i,:)')*0.01;
end

%% error of the Initial controller
Q = 10*eye(6);
R = 1*eye(3);
Cost_com = 0;
err_com = zeros(length(State_tra_com),6);
for i=1:length(State_tra)
    err_com(i,:) = get_error(State_tra_com(i,:))';
    ue = get_ue0(State_tra_com(i,:));
    Cost_com = Cost_com + (ue'*R*ue + err_com(i,:)*Q*err_com(i,:)')*0.01;
end

%% trajectory plot
fig1 = figure(1);
plot(State_tra(:,1), State_tra(:,2),'LineWidth',2)
hold on
plot(State_tra(:,6), State_tra(:,7),'--','LineWidth',2)
hold on
plot(State_tra_com(:,6), State_tra_com(:,7),':','LineWidth',2)
xlabel({'$x_L/x$ (m)'},'Interpreter','latex','FontSize',20)
ylabel({'$y_L/y$ (m)'},'Interpreter','latex','FontSize',20)
legend('Leader', 'ADP', 'Initial', 'Interpreter','latex','FontSize',20)
set(gca,'FontSize',15)   
H1 = axes('Position',[0.4,0.4,0.28,0.3]); 
plot(State_tra(1:200,1), State_tra(1:200,2),'LineWidth',2)
hold on
plot(State_tra(1:200,6), State_tra(1:200,7),'--','LineWidth',2)
hold on
plot(State_tra_com(1:200,6), State_tra_com(1:200,7),':','LineWidth',2)
saveas(fig1,'Plots/Tra_circle.eps', 'epsc');

fig2 = figure(2);
plot(T, err(:,1),'LineWidth',2)
hold on
plot(T, err_com(:,1),'LineWidth',1)
xlabel({'Time (s)'},'Interpreter','latex','FontSize',20)
ylabel({'$z_{1}$ (m)'},'Interpreter','latex','FontSize',20)
legend('ADP', 'Initial', 'Interpreter','latex','FontSize',20)
set(gca,'FontSize',15) 
% ylim([-15,5])
H1 = axes('Position',[0.25,0.6,0.28,0.3]); 
plot(T(1:100),  err(1:100,1),'LineWidth',2) 
hold on
plot(T(1:100),  err_com(1:100,1),'LineWidth',1) 
xlim([-0.1,1])
saveas(fig2,'Plots/z1_circle.eps', 'epsc');

fig3 = figure(3);
plot(T, err(:,2),'LineWidth',2)
hold on
plot(T, err_com(:,2),'LineWidth',1)
xlabel({'Time (s)'},'Interpreter','latex','FontSize',20)
ylabel({'$z_{2}$ (m)'},'Interpreter','latex','FontSize',20)
legend('ADP', 'Initial', 'Interpreter','latex','FontSize',20)
set(gca,'FontSize',15) 
H1 = axes('Position',[0.25,0.6,0.28,0.3]); 
plot(T(1:100), err(1:100,2),'LineWidth',2)
hold on
plot(T(1:100),  err_com(1:100,2),'LineWidth',1) 
xlim([-0.1,1])
saveas(fig3,'Plots/z2_circle.eps', 'epsc');

fig4 = figure(4);
plot(T, err(:,3),'LineWidth',2)
hold on
plot(T, err_com(:,3),'LineWidth',1)
xlabel({'Time (s)'},'Interpreter','latex','FontSize',20)
ylabel({'$e_{3}$ (rad)'},'Interpreter','latex','FontSize',20)
legend('ADP', 'Initial', 'Interpreter','latex','FontSize',20)
set(gca,'FontSize',15) 
H1 = axes('Position',[0.25,0.25,0.28,0.3]); 
plot(T(1:100), err(1:100,3),'LineWidth',2)
hold on
plot(T(1:100), err_com(1:100,3),'LineWidth',1)  
xlim([-0.1,1])
saveas(fig4,'Plots/e3_circle.eps', 'epsc');

% 
% % H1 = axes('Position',[0.2,0.22,0.28,0.3]); 
% % plot(T(1:100), e1(1:100),'LineWidth',2)
% % hold on
% % plot(T(1:100), e2(1:100),'LineWidth',2)
% % hold on
% % plot(T(1:100), e3(1:100),'LineWidth',2)              
% % ylim([-15,5]); 
% 
% figure(2)
% plot(T, e4,'LineWidth',2)
% hold on
% plot(T, e5,'LineWidth',2)
% hold on
% plot(T, e6,'LineWidth',2)
% xlabel({'Time (s)'},'Interpreter','latex','FontSize',23)
% ylabel({'$e_{4},e_{5},e_{6}$ (m/s or rad/s)'},'Interpreter','latex','FontSize',23)
% legend('$e_4$', '$e_5$', '$e_6$', 'Interpreter','latex','FontSize',23)
% set(gca,'FontSize',18)
% ylim([-10,5])

% H2 = axes('Position',[0.3,0.22,0.35,0.3]); 
% plot(T(1:100), e4(1:100),'LineWidth',2)
% hold on
% plot(T(1:100), e5(1:100),'LineWidth',2)
% hold on
% plot(T(1:100), e6(1:100),'LineWidth',2)              
% ylim([-10,5]); 

% e1 = X_test_Com(:,1) + s*sin(X_test_Com(:,3)) - X_test_Com(:,6) - d*cos(X_test_Com(:,8));
% e2 = X_test_Com(:,2) - s*cos(X_test_Com(:,3)) - X_test_Com(:,7) - d*sin(X_test_Com(:,8));
% e3 = X_test_Com(:,3) - (X_test_Com(:,8) + gamma);
% e4 = X_test_Com(:,4) - X_test_Com(:,9);
% e5 = -X_test_Com(:,10);
% e6 = X_test_Com(:,5) - X_test_Com(:,11);
% z1 =  cos(X_test_Com(:,8)+gamma).*e1 + sin(X_test_Com(:,8)+gamma).*e2;
% z2 = -sin(X_test_Com(:,8)+gamma).*e1 + cos(X_test_Com(:,8)+gamma).*e2;
% 
% figure(3)
% plot(T_test, e1,'LineWidth',2)
% hold on
% plot(T_test, e2,'LineWidth',2)
% hold on
% plot(T_test, e3,'LineWidth',2)
% 
% xlabel({'Time (s)'},'Interpreter','latex','FontSize',23)
% ylabel({'$e_{1},e_{2},e_{3}$ (m or rad)'},'Interpreter','latex','FontSize',23)
% legend('$e_1$', '$e_2$', '$e_3$', 'Interpreter','latex','FontSize',23)
% set(gca,'FontSize',18)   
% 
% H3 = axes('Position',[0.2,0.22,0.28,0.3]); 
% plot(T_test(1:100), e1(1:100),'LineWidth',2)
% hold on
% plot(T_test(1:100), e2(1:100),'LineWidth',2)
% hold on
% plot(T_test(1:100), e3(1:100),'LineWidth',2)              
% ylim([-15,5]); 
% 
% figure(4)
% plot(T_test, e4,'LineWidth',2)
% hold on
% plot(T_test, e5,'LineWidth',2)
% hold on
% plot(T_test, e6,'LineWidth',2)
% xlabel({'Time (s)'},'Interpreter','latex','FontSize',23)
% ylabel({'$e_{4},e_{5},e_{6}$ (m/s or rad/s)'},'Interpreter','latex','FontSize',23)
% legend('$e_4$', '$e_5$', '$e_6$', 'Interpreter','latex','FontSize',23)
% set(gca,'FontSize',18)
% 
% H4 = axes('Position',[0.3,0.22,0.35,0.3]); 
% plot(T_test(1:100), e4(1:100),'LineWidth',2)
% hold on
% plot(T_test(1:100), e5(1:100),'LineWidth',2)
% hold on
% plot(T_test(1:100), e6(1:100),'LineWidth',2)              
% ylim([-10,5]); 
% 
% figure(1)
% plot(t_save, Data_save(:,1:5), 'LineWidth',2)
% 
% xlabel({'Time (s)'},'Interpreter','latex','FontSize',23)
% ylabel({'Leading Vehicle'},'Interpreter','latex','FontSize',23)
% legend('$x_L$', '$y_L$', '$\phi_L$', '$V_L$','$\omega_L$', 'Interpreter','latex')
% 
% figure(2)
% plot(t_save, Data_save(:,6:11), 'LineWidth',2)
% 
% xlabel({'Time (s)'},'Interpreter','latex','FontSize',23)
% ylabel({'Following Vehicle'},'Interpreter','latex','FontSize',23)
% legend('$x$', '$y$', '$\phi$', '$V_x$', '$V_y$', '$\omega_L$', 'Interpreter','latex')



% % 
% % figure(3)
% % plot(1:1:12, J(1,:),'-gs','LineWidth',2, 'MarkerSize',10,'MarkerEdgeColor','b','MarkerFaceColor',[0.5,0.5,0.5])
% % xlabel({'Iteration'},'Interpreter','latex','FontSize',23)
% % ylabel({'$J$'},'Interpreter','latex','FontSize',23)
% % set(gca,'FontSize',18)