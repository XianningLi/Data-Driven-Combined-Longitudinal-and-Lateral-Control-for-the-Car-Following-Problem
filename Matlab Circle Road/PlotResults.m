%% Plot the Training Results
NORM_SU = [];
for i=2:iter
    NORM_SU = [NORM_SU, norm(SU(:,:,i))];
end
figure(1)
% plot(1:length(NORM_SU), NORM_SU,'-o','LineWidth',2,'MarkerSize',6,'MarkerEdgeColor','b','MarkerFaceColor',[0.5,0.5,0.5])
plot(1:length(NORM_SU), NORM_SU,'-o','LineWidth',2)
xlabel({'Number of Iteration'},'Interpreter','latex','FontSize',18)
legend('norm of weights of control policy', 'Interpreter','latex','FontSize',15)
set(gca,'FontSize',15)
xlim([1 length(NORM_SU)])
% xticks([1:5:50,53])

NORM_SV = [];
for i=1:iter-1
    NORM_SV = [NORM_SV, norm(SV(:,i))];
end
figure(2)
% plot(1:length(NORM_SV), NORM_SV,'-o','LineWidth',2,'MarkerSize',6,'MarkerEdgeColor','b','MarkerFaceColor',[0.5,0.5,0.5])
plot(1:length(NORM_SV), NORM_SV,'-o','LineWidth',2)
xlabel({'Number of Iteration'},'Interpreter','latex','FontSize',18)
legend('norm of weights of value function', 'Interpreter','latex','FontSize',15)
set(gca,'FontSize',15)
xlim([1 length(NORM_SV)])
% xticks([1:5:50,53])

figure(3)
plot(x_save(:,1), x_save(:,2),'LineWidth',2)
hold on
plot(x_save(:,6), x_save(:,7),'LineWidth',2)
legend('Leading Vehicle', 'Following Vehicle', 'Interpreter','latex','FontSize',15)
xlabel({'$x_L/x (m)$'},'Interpreter','latex','FontSize',18)
ylabel({'$y_L/y (m)$'},'Interpreter','latex','FontSize',18)
set(gca,'FontSize',15)

figure(4)
plot(t_save, x_save(:,1),'LineWidth',2)
hold on
plot(t_save, x_save(:,6),'LineWidth',2)
legend('Leading Vehicle', 'Following Vehicle', 'Interpreter','latex','FontSize',15)
xlabel({'Time(s)'},'Interpreter','latex','FontSize',18)
ylabel({'$x_L/x(m)$'},'Interpreter','latex','FontSize',18)
set(gca,'FontSize',15)

figure(5)
plot(t_save, x_save(:,2),'LineWidth',2)
hold on
plot(t_save, x_save(:,7),'LineWidth',2)
legend('Leading Vehicle', 'Following Vehicle', 'Interpreter','latex','FontSize',15)
xlabel({'Time(s)'},'Interpreter','latex','FontSize',18)
ylabel({'$y_L/y(m)$'},'Interpreter','latex','FontSize',18)
set(gca,'FontSize',15)

figure(6)
plot(t_save, x_save(:,3),'LineWidth',2)
hold on
plot(t_save, x_save(:,8),'LineWidth',2)
legend('Leading Vehicle', 'Following Vehicle', 'Interpreter','latex','FontSize',15)
xlabel({'Time(s)'},'Interpreter','latex','FontSize',18)
ylabel({'$\varphi_L/\varphi(rad)$'},'Interpreter','latex','FontSize',18)
set(gca,'FontSize',15)

% radius =  100;
% kappa = 1/radius;
% d = 20;
% gamma = atan(kappa*d);
% 
% if kappa == 0
%     s = 0;
% else
%     s = (-1 + sqrt(1+kappa^2*d^2))/kappa;
% end
% 
% e1 = x_save(:,1) + s*sin(x_save(:,3)) - x_save(:,6) - d*cos(x_save(:,8));
% e2 = x_save(:,2) - s*cos(x_save(:,3)) - X_test(:,7) - d*sin(X_test(:,8));
% e3 = x_save(:,3) - (x_save(:,8) + gamma);



% radius =  100;
% kappa = 1/radius;
% d = 20;
% gamma = atan(kappa*d);
% 
% if kappa == 0
%     s = 0;
% else
%     s = (-1 + sqrt(1+kappa^2*d^2))/kappa;
% end
% 
% e1 = X_test(:,1) + s*sin(X_test(:,3)) - X_test(:,6) - d*cos(X_test(:,8));
% e2 = X_test(:,2) - s*cos(X_test(:,3)) - X_test(:,7) - d*sin(X_test(:,8));
% e3 = X_test(:,3) - (X_test(:,8) + gamma);
% e4 = X_test(:,4) - X_test(:,9);
% e5 = -X_test(:,10);
% e6 = X_test(:,5) - X_test(:,11);
% z1 =  cos(X_test(:,8)+gamma).*e1 + sin(X_test(:,8)+gamma).*e2;
% z2 = -sin(X_test(:,8)+gamma).*e1 + cos(X_test(:,8)+gamma).*e2;
% 
% figure(1)
% plot(T_test, e1,'LineWidth',2)
% hold on
% plot(T_test, e2,'LineWidth',2)
% hold on
% plot(T_test, e3,'LineWidth',2)
% ylim([-15,5])
% 
% xlabel({'Time (s)'},'Interpreter','latex','FontSize',23)
% ylabel({'$e_{1},e_{2},e_{3}$ (m or rad)'},'Interpreter','latex','FontSize',23)
% legend('$e_1$', '$e_2$', '$e_3$', 'Interpreter','latex','FontSize',23)
% set(gca,'FontSize',18)   
% 
% H1 = axes('Position',[0.2,0.22,0.28,0.3]); 
% plot(T_test(1:100), e1(1:100),'LineWidth',2)
% hold on
% plot(T_test(1:100), e2(1:100),'LineWidth',2)
% hold on
% plot(T_test(1:100), e3(1:100),'LineWidth',2)              
% ylim([-15,5]); 
% 
% figure(2)
% plot(T_test, e4,'LineWidth',2)
% hold on
% plot(T_test, e5,'LineWidth',2)
% hold on
% plot(T_test, e6,'LineWidth',2)
% xlabel({'Time (s)'},'Interpreter','latex','FontSize',23)
% ylabel({'$e_{4},e_{5},e_{6}$ (m/s or rad/s)'},'Interpreter','latex','FontSize',23)
% legend('$e_4$', '$e_5$', '$e_6$', 'Interpreter','latex','FontSize',23)
% set(gca,'FontSize',18)
% ylim([-10,5])
% 
% H2 = axes('Position',[0.3,0.22,0.35,0.3]); 
% plot(T_test(1:100), e4(1:100),'LineWidth',2)
% hold on
% plot(T_test(1:100), e5(1:100),'LineWidth',2)
% hold on
% plot(T_test(1:100), e6(1:100),'LineWidth',2)              
% ylim([-10,5]); 
% 
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
% 
% 
% % 
% % figure(3)
% % plot(1:1:12, J(1,:),'-gs','LineWidth',2, 'MarkerSize',10,'MarkerEdgeColor','b','MarkerFaceColor',[0.5,0.5,0.5])
% % xlabel({'Iteration'},'Interpreter','latex','FontSize',23)
% % ylabel({'$J$'},'Interpreter','latex','FontSize',23)
% % set(gca,'FontSize',18)