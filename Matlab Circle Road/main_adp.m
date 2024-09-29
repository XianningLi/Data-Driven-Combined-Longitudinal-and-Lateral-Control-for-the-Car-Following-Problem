clear all

%% parameters of the road
radius =  100;
kappa = 1/radius;
ts = 1;
ds = 5;
d = ts*20 + ds;
gamma = atan(kappa*d);
if kappa == 0
    s = 0;
else
    s = (-1 + sqrt(1+kappa^2*d^2))/kappa;
end

%% Initial state and parameters of the system
X0 = [radius, 0, pi/2, 20, 20*kappa,  radius*cos(gamma)+5, -radius*sin(gamma)+5, pi/2 - gamma+0.2, 20-5, 0, 20*kappa];

un = 3;
xn = 11;
Num_bas_val = length(get_basis_Val(X0));
Num_bas_ue = length(get_basis_ue(X0));
Num_basisfV = length(basis_f_V(X0));

%% parameters of the learning algorithm
Q = 10*eye(6);
R = 1*eye(un);

T = 0.01;
expl_noise_freq = (rand(un,500)-.5)*500; 
Su = 0;
Sf = 0;
g = 0;

%% Data Collection for phase-one
x_save = [];
t_save = [];

ifLearned = 0;
ifPhase1 = 1;
step = 0;
N1 = 50;
X_aug = [X0, zeros(1,Num_basisfV*3), kron(zeros(1,un),X0(9:11))];

Dxx4=[];
Ixx4=[];
Ixu4=[];
Dxx5=[];
Ixx5=[];
Ixu5=[];
Dxx6=[];
Ixx6=[];
Ixu6=[];

for i=1:N1
    % Simulation the system and at the same time collect online info.
    [t,X_aug] = ode45(@(t,X_aug)aug_sys(t, X_aug, ifPhase1, ifLearned, expl_noise_freq, Q, R, Su, Sf, g, step), ...
        [(i-1)*T,i*T],X_aug(end,:));
    
    %Append new data to the data matrices
    Dxx4=[Dxx4;0.5*(X_aug(end,9)^2 - X_aug(1,9)^2)];
    Ixx4=[Ixx4;X_aug(end,xn+1:xn+Num_basisfV)-X_aug(1,xn+1:xn+Num_basisfV)];
    Ixu4=[Ixu4;X_aug(end,xn+Num_basisfV+1:xn+Num_basisfV+un)-X_aug(1,xn+Num_basisfV+1:xn+Num_basisfV+un)];
    
    Dxx5=[Dxx5;0.5*(X_aug(end,10)^2 - X_aug(1,10)^2)];
    Ixx5=[Ixx5;X_aug(end,xn+(un+Num_basisfV)+1:xn+(un+Num_basisfV)+Num_basisfV)-X_aug(1,xn+(un+Num_basisfV)+1:xn+(un+Num_basisfV)+Num_basisfV)];
    Ixu5=[Ixu5;X_aug(end,xn+(un+Num_basisfV)+Num_basisfV+1:xn+(un+Num_basisfV)+Num_basisfV+un)-X_aug(1,xn+(un+Num_basisfV)+Num_basisfV+1:xn+(un+Num_basisfV)+Num_basisfV+un)];
    
    Dxx6=[Dxx6;0.5*(X_aug(end,11)^2 - X_aug(1,11)^2)];
    Ixx6=[Ixx6;X_aug(end,xn+(un+Num_basisfV)+(un+Num_basisfV)+1:xn+(un+Num_basisfV)+(un+Num_basisfV)+Num_basisfV)-X_aug(1,xn+(un+Num_basisfV)+(un+Num_basisfV)+1:xn+(un+Num_basisfV)+(un+Num_basisfV)+Num_basisfV)];
    Ixu6=[Ixu6;X_aug(end,xn+(un+Num_basisfV)+(un+Num_basisfV)+Num_basisfV+1:xn+(un+Num_basisfV)+(un+Num_basisfV)+Num_basisfV+un)-X_aug(1,xn+(un+Num_basisfV)+(un+Num_basisfV)+Num_basisfV+1:xn+(un+Num_basisfV)+(un+Num_basisfV)+Num_basisfV+un)];
    
    % Keep track of the system trajectories
    x_save=[x_save;X_aug(end,1:xn)];
    t_save=[t_save;t(end)]; 
end

%% Learning for Phase-one
param1 = pinv([Ixx4,Ixu4])*Dxx4;
param2 = pinv([Ixx5,Ixu5])*Dxx5;
param3 = pinv([Ixx6,Ixu6])*Dxx6;

Sf = [param1(1:Num_basisfV)';param2(1:Num_basisfV)';param3(1:Num_basisfV)'];
g = [param1(Num_basisfV+1:end)';param2(Num_basisfV+1:end)';param3(Num_basisfV+1:end)'];

%% Data Collection for phase-two
% The initial controller cannot be expressed as the linear
% combiniation of the given basis function. 

N2 = 1000;
ifLearned = 0;
ifPhase1 = 0;
step = 0;

X_aug0 = [x_save(end,:), zeros(1,Num_bas_ue*un), 0];
Data0 = X_aug0;
for i=N1+1:N1+N2
    [t,X_aug] = ode45(@(t,X_aug)aug_sys(t, X_aug, ifPhase1, ifLearned, expl_noise_freq, Q, R, Su, Sf, g, step), ...
        [i*T: 0.001: (i+1)*T],X_aug0); 
    X = X_aug(end,1:xn);
    int_nu_bu = X_aug(end,xn+1:xn+Num_bas_ue*un);
    int_cost = X_aug(end, xn+Num_bas_ue*un+1);
    Data0 = [Data0; [X, int_nu_bu, int_cost]];
    X_aug0 = [X, zeros(1,Num_bas_ue*un), 0];
end

step = 2;
X_aug0 = [x_save(end,:), zeros(1,Num_bas_ue^2), zeros(1,Num_bas_ue*un), 0];
Data_save = X_aug0;

for i=N1+1:N1+N2
    [t,X_aug] = ode45(@(t,X_aug)aug_sys(t, X_aug, ifPhase1, ifLearned, expl_noise_freq, Q, R, Su, Sf, g, step), ...
        [i*T: 0.001: (i+1)*T],X_aug0);    
    
    X = X_aug(end,1:xn);
    int_bu_bu = X_aug(end,xn+1:xn+Num_bas_ue^2);
    int_ue_bu = X_aug(end,xn+Num_bas_ue^2+1:xn+Num_bas_ue^2+un*Num_bas_ue);
    int_Qe = X_aug(end,xn+Num_bas_ue^2+un*Num_bas_ue+1);
    Data_save = [Data_save; [X, int_bu_bu, int_ue_bu, int_Qe]];
    X_aug0 = [X, zeros(1,Num_bas_ue^2), zeros(1,Num_bas_ue*un), 0];
    
    x_save=[x_save;X_aug(end,1:xn)];
    t_save=[t_save; t(end)];
end

%% Learning process for phase-two
%Start from the initial controller
Num_iter = 100;
SU = zeros(un, Num_bas_ue, Num_iter);
SV = zeros(Num_bas_val, Num_iter);

Diff_base_Val = zeros(N2,Num_bas_val);
parfor i=1:N2
    Diff_base_Val(i,:) = get_basis_Val(Data0(i+1,1:xn))'-get_basis_Val(Data0(i,1:xn))';
end
Int_nu_bu = Data0(2:end,xn+1:xn+Num_bas_ue*un);
Int_cost = Data0(2:end,xn+Num_bas_ue*un+1);
[ZZ,RankFlag] = chol([Diff_base_Val, 2*Int_nu_bu]'*[Diff_base_Val, 2*Int_nu_bu]);
RankFlag
prmt = -pinv([Diff_base_Val, 2*Int_nu_bu])*Int_cost;
Sv = prmt(1:Num_bas_val);
Su = reshape(prmt(Num_bas_val+1:end),[Num_bas_ue,un])';
SU(:,:,2) = Su;
Sv(:,:,1) = Sv;

Su_pre = zeros(un, Num_bas_ue);
iter = 2;
while norm(Su-Su_pre)>10^(-3) && iter<Num_iter
    iter = iter + 1;
    Diff_base_Val = zeros(N2, Num_bas_val);
    Int_nu_bu = zeros(N2, Num_bas_ue*un);
    Int_cost = zeros(N2,1);
    parfor i=1:N2
        Diff_base_Val(i,:) = get_basis_Val(Data_save(i+1,1:xn))'-get_basis_Val(Data_save(i,1:xn))';
        int_bu_bu = Data_save(i+1,xn+1:xn+Num_bas_ue^2);
        int_ue_bu = Data_save(i+1,xn+Num_bas_ue^2+1:xn+Num_bas_ue^2+un*Num_bas_ue);
        Qe = Data_save(i+1,xn+Num_bas_ue^2+un*Num_bas_ue+1);
        Ru = int_bu_bu*reshape(Su'*R*Su,[Num_bas_ue^2,1]);
        Int_cost(i,:) = Qe+Ru;
        Int_nu_bu(i,:) = int_ue_bu-int_bu_bu*kron(Su'*R',eye(Num_bas_ue));
    end 
    [ZZ,RankFlag] = chol([Diff_base_Val, 2*Int_nu_bu]'*[Diff_base_Val, 2*Int_nu_bu]);
    RankFlag
    prmt = -pinv([Diff_base_Val, 2*Int_nu_bu])*Int_cost;
    Sv = prmt(1:Num_bas_val);
    Su_pre = Su;
    Su = reshape(prmt(Num_bas_val+1:end),[Num_bas_ue,un])';
    SU(:,:,iter) = Su;
    SV(:,iter-1) = Sv;
end

