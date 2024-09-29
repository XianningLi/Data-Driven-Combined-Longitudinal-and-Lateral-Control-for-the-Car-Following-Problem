function ue0 = get_ue0(X) 
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

K1 = 1*eye(3);
K2 = 4*eye(3);

%% Get the state of the vehicles
x_L =  X(1);
y_L =  X(2);
phi_L =  X(3); 
V_L =  X(4); 
omega_L =  X(5);

x =  X(6);
y =  X(7);
phi =  X(8); 
V_x =  X(9); 
V_y =  X(10);
omega =  X(11);

%% parameters of the path
radius =  100;
kappa = 1/radius;
ts = 1;
ds = 5;
d = ts*V_x + ds;
gamma = atan(kappa*d);

if kappa == 0
    s = 0;
else
    s = (-1 + sqrt(1+kappa^2*d^2))/kappa;
end

%% Get the error state
beta = atan(V_y/V_x);
e1 = x_L + s*sin(phi_L) - x - d*cos(phi);
e2 = y_L - s*cos(phi_L) - y - d*sin(phi);
e3 = phi_L - (phi + gamma);
e4 = V_L - V_x;
e5 = -V_y;
e6 = omega_L - omega;

e123 = [e1;e2;e3];
e456 = [e4;e5;e6];
%% Calculate Controller
f_ep = [V_L*cos(phi_L) + s*cos(phi_L)*omega_L - V_L*cos(phi_L - gamma - e3) + d*sin(phi_L - gamma - e3)*omega_L;...
        V_L*sin(phi_L) + s*sin(phi_L)*omega_L - V_L*sin(phi_L - gamma - e3) - d*cos(phi_L - gamma - e3)*omega_L;...
        0];
f_ep_dot = [-V_L*sin(phi_L)*omega_L - s*sin(phi_L)*omega_L^2 + V_L*sin(phi_L-gamma-e3)*(omega_L - e6) + d*cos(phi_L-gamma-e3)*omega_L*(omega_L - e6);...
             V_L*cos(phi_L)*omega_L + s*cos(phi_L)*omega_L^2 - V_L*cos(phi_L-gamma-e3)*(omega_L - e6) + d*sin(phi_L-gamma-e3)*omega_L*(omega_L - e6);...
             0];
g_ep = [ cos(phi_L - gamma - e3), -sin(phi_L - gamma - e3), -d*sin(phi_L - gamma - e3);...
         sin(phi_L - gamma - e3),  cos(phi_L - gamma - e3),  d*cos(phi_L - gamma - e3);...
           0,                       0,                       1];
g_ep_inv = [  cos(phi_L-gamma-e3), sin(phi_L-gamma-e3),  0;...
             -sin(phi_L-gamma-e3), cos(phi_L-gamma-e3), -d;...
             0,                       0,                       1];
g_ep_inv_dot = [-sin(phi_L-gamma-e3),  cos(phi_L-gamma-e3),   0;...
                -cos(phi_L-gamma-e3),  -sin(phi_L-gamma-e3),  0;...
                0,                       0,                   0]*(omega_L - e6);
e123_dot = f_ep + g_ep*e456;
alphae = g_ep_inv*( -f_ep - K1*e123);
alphae_dot = g_ep_inv_dot*( -f_ep - K1*[e1;e2;e3]) + g_ep_inv*( -f_ep_dot - K1*e123_dot);

% f_ev = [-V_y*omega + Ca/M*V_x^2 - Ca/M*V_L^2;...
%          V_x*omega - V_L*omega_L + (Cr + Cf)/M*V_y/V_x + (Cr*lr - Cf*lf)/M*(omega_L/V_L - omega/V_x);...
%          -(Cr*lr - Cf*lf)/Iz*V_y/V_x + (Cf*lf^2 + Cr*lr^2)/Iz*(omega/V_x - omega_L/V_L)];
f_ev = [-V_y*omega + Ca/M*V_x^2 - Ca/M*V_L^2;...
         V_x*omega - V_L*omega_L + (Cr + Cf)/M*V_y/V_x + (Cr*lr - Cf*lf)/M*(omega_L/V_L - omega/V_x);...
         -(Cr*lr - Cf*lf)/Iz*V_y/V_x + (Cf*lf^2 + Cr*lr^2)/Iz*(omega/V_x- omega_L/V_L)];
          
g_v = [2*k0/M,       0,         2*k0/M;...
       0,            Cf/M,      0;...
       - ls*2*k0/Iz, Cf*lf/Iz,  ls*2*k0/Iz;];


ue0 = inv(g_v)*(f_ev - alphae_dot + K2*(e456 - alphae) + g_ep'*e123);
end