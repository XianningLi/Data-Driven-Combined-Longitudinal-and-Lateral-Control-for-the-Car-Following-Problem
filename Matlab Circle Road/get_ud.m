function ud = get_ud(X) 
%% Parameters of the following vehicle
M = 1360;
Iz = 1993;
lf = 1.45;
lr = 1.06;
ls = 0.71;
Ca = 0.5;
k0 = 460/0.33;
Cf = 1.51*10^5;
Cr = 1.46*10^5;


%% Get the state of the vehicles
x_L =  X(1);
y_L =  X(2);
phi_L =  X(3); 
V_L =  X(4); 
omega_L =  X(5);
     
g_v = [2*k0/M,       0,         2*k0/M;...
       0,            Cf/M,      0;...
       - ls*2*k0/Iz, Cf*lf/Iz,  ls*2*k0/Iz;];

ud = inv(g_v)*[Ca/M*V_L^2;...
               V_L*omega_L - (Cr*lr - Cf*lf)/M*omega_L/V_L;...
               (Cf*lf^2 + Cr*lr^2)/Iz*omega_L/V_L];

end