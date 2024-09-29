
kappa = 1/3;
d = 1;
gamma = atan(kappa*d);
if kappa == 0
s = 0;
else
s = (-1 + sqrt(1+kappa^2*d^2))/kappa;
end

X0 = [3,0,pi/2,24,8,  3*cos(gamma)+ 1, -3*sin(gamma)+1,pi/2 - gamma,24,0, 8, 0]';
% X0 = [3,0,pi/2,3,1,  10,20,pi/2,10,2,10]';

un = 3;
Num_X = 11;
Num_bas_ue = length(get_basis_ue(X0));
ifLearned = 1;
ifPhase1 = 0;
expl_noise_freq = (rand(un,200)-.5)*20;
step = 1;
X_aug = X0';
X_aug(6:11) = X_aug(6:11) + (rand(1,Num_X - 5)-.5)*20;

x_save = [];
t_save = [];
T = 0.01;
N = 10000;



[t_save,x_save] = ode45(@(t,X_aug)aug_sys(t, X_aug, ifPhase1, ifLearned, expl_noise_freq, Su, step), ...
    0: 0.001: 20, X0);

