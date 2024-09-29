radius =  100;
kappa = 1/radius;
d = 20;
gamma = atan(kappa*d);
if kappa == 0
s = 0;
else
s = (-1 + sqrt(1+kappa^2*d^2))/kappa;
end

%% Initial state and parameters of the system
X0 = [radius, 0, pi/2, 10, 10*kappa,  radius*cos(gamma)+3, -radius*sin(gamma)+3, pi/2 - gamma, 5, 0, 5*kappa];
X0 = [X0,0];

ifPhase1 = 0;
ifLearned = 1;
expl_noise_freq = 0;
step = 0;

[T_test,X_test] = ode45(@(t,X_aug)aug_sys(t, X_aug, ifPhase1, ifLearned, expl_noise_freq, Q, R, Su, Sf, g, step), ...
    0: 0.01: 50, X0);


Err = zeros(size(X_test,1), 6);
for i = 1:size(X_test,1)
   Err(i,:) = get_error(X_test(i,1:11)); 
end

plot(T_test, Err)


