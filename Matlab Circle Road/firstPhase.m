clear all
%% parameters of the road
radius =  100;
kappa = 1/radius;
d = 20;
gamma = atan(kappa*d);
if kappa == 0
    s = 0;
else
    s = (-1 + sqrt(1+kappa^2*d^2))/kappa;
end
%% 
un = 3;
Num_X = 11;
Num_V = 3;
Num_basisfV = 5;

X0 = [radius, 0, pi/2, 20, 20*kappa,  radius*cos(gamma)+5, -radius*sin(gamma)+5, pi/2 - gamma+0.2, 20-5, 0, 20*kappa];
X_aug = [X0; kron(zeros(1,Num_basisfV),X0(9:11)')'; kron(zeros(1,un),X0(9:11)')']';

T  = 0.01; 
N = 50;

Dxx4=[];
Ixx4=[];
Ixu4=[];
Dxx5=[];
Ixx5=[];
Ixu5=[];
Dxx6=[];
Ixx6=[];
Ixu6=[];

ifLearned = 0;
ifPhase1 = 1;
expl_noise_freq = (rand(un,500)-.5)*500; 
Su = 0;
step = 0;

Q = 10*eye(6);
R = 1*eye(un);

x_save=[];
t_save=[];

%%
for i=1:N
    % Simulation the system and at the same time collect online info.
    [t,X_aug] = ode45(@(t,X_aug)aug_sys(t, X_aug, ifPhase1, ifLearned, expl_noise_freq, Q, R, Su, step), ...
        [(i-1)*T,i*T],X_aug(end,:));
    
    %Append new data to the data matrices
    Dxx4=[Dxx4;0.5*(X_aug(end,9)^2 - X_aug(1,9)^2)];
    Ixx4=[Ixx4;X_aug(end,Num_X+1:Num_X+Num_basisfV)-X_aug(1,Num_X+1:Num_X+Num_basisfV)];
    Ixu4=[Ixu4;X_aug(end,Num_X+Num_basisfV+1:Num_X+Num_basisfV+un)-X_aug(1,Num_X+Num_basisfV+1:Num_X+Num_basisfV+un)];
    
    Dxx5=[Dxx5;0.5*(X_aug(end,10)^2 - X_aug(1,10)^2)];
    Ixx5=[Ixx5;X_aug(end,Num_X+(un+Num_basisfV)+1:Num_X+(un+Num_basisfV)+Num_basisfV)-X_aug(1,Num_X+(un+Num_basisfV)+1:Num_X+(un+Num_basisfV)+Num_basisfV)];
    Ixu5=[Ixu5;X_aug(end,Num_X+(un+Num_basisfV)+Num_basisfV+1:Num_X+(un+Num_basisfV)+Num_basisfV+un)-X_aug(1,Num_X+(un+Num_basisfV)+Num_basisfV+1:Num_X+(un+Num_basisfV)+Num_basisfV+un)];
    
    Dxx6=[Dxx6;0.5*(X_aug(end,11)^2 - X_aug(1,11)^2)];
    Ixx6=[Ixx6;X_aug(end,Num_X+2*(un+Num_basisfV)+1:Num_X+2*(un+Num_basisfV)+Num_basisfV)-X_aug(1,Num_X+2*(un+Num_basisfV)+1:Num_X+2*(un+Num_basisfV)+Num_basisfV)];
    Ixu6=[Ixu6;X_aug(end,Num_X+2*(un+Num_basisfV)+Num_basisfV+1:Num_X+2*(un+Num_basisfV)+Num_basisfV+un)-X_aug(1,Num_X+2*(un+Num_basisfV)+Num_basisfV+1:Num_X+2*(un+Num_basisfV)+Num_basisfV+un)];
    
    % Keep track of the system trajectories
    x_save=[x_save;X_aug];
    t_save=[t_save;t];
    
end

Lambda = [Ixx4, Ixx5, Ixx6, Ixu4, Ixu5, Ixu6];
XI = Dxx4+Dxx5+Dxx6;

param = pinv(Lambda)*XI;
plot(t_save,x_save(:,1:11))


