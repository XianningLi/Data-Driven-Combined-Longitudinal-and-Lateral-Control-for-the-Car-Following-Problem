clear all

un = 3;
Num_X = 11;
Num_V = 3;
Num_basisfV = 5;
% expl_noise_freq = (rand(un,5)-.5)*10;
expl_noise_freq = [1;2;3];

X0 = [3,0,pi/2,3,1,  3,0,pi/2,3,0.3,1]';
X_aug = [X0; kron(zeros(1,Num_basisfV),X0(9:11)')'; kron(zeros(1,un),X0(9:11)')']';

T  = 0.01; 
N = 1000;

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

x_save=[];
t_save=[];

for i=1:N
    % Simulation the system and at the same time collect online info.
    [t,X_aug] = ode45(@(t,X_aug)aug_sys(t,X_aug,ifLearned,expl_noise_freq), ...
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


