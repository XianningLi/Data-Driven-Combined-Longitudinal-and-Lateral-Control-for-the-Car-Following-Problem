function dX_aug = test_aug_sys(t, X_aug, u)

Q = 10*eye(6);
R = 0.1*eye(3);

X = X_aug(1:11);
e = get_error(X);
% u = get_ue0(X) + get_ud(X);


dX_aug = actual_sys(X,u);
end

function dX = actual_sys(X,u)
%DINAMICS the dynamic model of the leading vehicle and the following
%vehicle
%   X(1) = x_L, X(2) = y_L, X(3) = \phi_L, X(4) = V_L, X(5) = \omega
%   X(6) = x, X(7) = y, X(8) = \phi, X(9) = V_x, X(10) = V_x, X(11) = V_y
%   
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
%% parameters of the path
kappa = 1/3;
d = 1;
gamma = atan(kappa*d);

if kappa == 0
    s = 0;
else
    s = (-1 + sqrt(1+kappa^2*d^2))/kappa;
end
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

%% Get the error state
beta = atan(V_y/V_x);
e1 = x_L + s*sin(phi_L) - x - d*cos(phi);
e2 = y_L - s*cos(phi_L) - y - d*sin(phi);
e3 = phi_L - (phi + gamma);
e4 = V_L - V_x;
e5 = -V_y;
e6 = omega_L - omega;

%% Dynamics of the leading vehicle
x_L_dot = V_L*cos(phi_L);
y_L_dot = V_L*sin(phi_L);
phi_L_dot = omega_L;
V_L_dot = 0;
omega_L_dot = 0;

%%
u1 = u(1);
u2 = u(2);
u3 = u(3);
%% Dynamics of the following vehicle
x_dot = V_x*cos(phi) - V_y*sin(phi);
y_dot = V_x*sin(phi) + V_y*cos(phi);
phi_dot = omega;
V_x_dot = V_y*omega - Ca/M*V_x^2 + 2*k0/M*u1 + 2*k0/M*u3;
V_y_dot = -V_x*omega - (Cf + Cr)/M*V_y/V_x + (Cr*lr - Cf*lf)/M*omega/V_x  + Cf/M*u2;
omega_dot = (Cr*lr - Cf*lf)/Iz*V_y/V_x - (Cf*lf^2 + Cr*lr^2)/Iz*omega/V_x - ls*2*k0/Iz*u1 + ls*2*k0/Iz*u3 + Cf*lf/Iz*u2;

dX = [x_L_dot; y_L_dot; phi_L_dot; V_L_dot; omega_L_dot; x_dot; y_dot; phi_dot; V_x_dot; V_y_dot; omega_dot];

end

