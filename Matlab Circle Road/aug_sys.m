function dX_aug = aug_sys(t, X_aug, ifPhase1, ifLearned, expl_noise_freq, Q, R, Su, Sf, g, step)

X = X_aug(1:11);
e = get_error(X);
% ud = get_ud(X);
V_L = X(4);
omega_L = X(5);

if ~ifLearned
    if ifPhase1
        u = 0.001*sum(sin(expl_noise_freq*t),2);
        basis_f = basis_f_V(X);
        dX = actual_sys(X,u);
        dXX4 = basis_f*X(9);
        duX4 = u*X(9);
        dXX5 = basis_f*X(10);
        duX5 = u*X(10);
        dXX6 = basis_f*X(11);
        duX6 = u*X(11);
        dX_regu = [dXX4;duX4;dXX5;duX5;dXX6;duX6];
        dX_aug = [dX;dX_regu];
    else
        ud = pinv(g)*(-Sf*basis_f_V([X(1:8);V_L;0;omega_L]));
        ue = 0.001*sum(sin(expl_noise_freq*t),2)+get_ue0(X); %04/02/2022
        u = ue + ud;
        basis_ue = get_basis_ue(X);
        if step < 1 
            uei = get_ue0(X);
            nui = ue - uei; 
            dX = actual_sys(X,u);
            d_nu_bu = kron(R*nui,basis_ue);
            dcost = e'*Q*e + uei'*R*uei;
            dX_aug = [dX; d_nu_bu; dcost];                      
        else
            dX = actual_sys(X,u);
            d_bu_bu = kron(basis_ue,basis_ue);
            d_ue_bu = kron(R*ue,basis_ue);
            dQe = e'*Q*e;
            dX_aug = [dX; d_bu_bu; d_ue_bu; dQe];
        end
    end
else
    ud = pinv(g)*(-Sf*basis_f_V([X(1:8);V_L;0;omega_L]));
    ue = Su*get_basis_ue(X);
%     ue = get_ue0(X);
%     ue = 0.002*sum(sin(expl_noise_freq*t),2);
    u = ue + ud;
    Qe = e'*Q*e;
    Ru = ue'*R*ue;
    dX = actual_sys(X,u);
    dcost = Qe + Ru;
    dX_aug = [dX;dcost];    
end

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
% V_x_dot = V_y*omega - Ca/M*V_x^2 + 2*k0/M*u1 + 2*k0/M*u3;
% V_y_dot = -V_x*omega - (Cf + Cr)/M*V_y/V_L + (Cr*lr - Cf*lf)/M*omega/V_L  + Cf/M*u2;
% omega_dot = (Cr*lr - Cf*lf)/Iz*V_y/V_L - (Cf*lf^2 + Cr*lr^2)/Iz*omega/V_L - ls*2*k0/Iz*u1 + ls*2*k0/Iz*u3 + Cf*lf/Iz*u2;

dX = [x_L_dot; y_L_dot; phi_L_dot; V_L_dot; omega_L_dot; x_dot; y_dot; phi_dot; V_x_dot; V_y_dot; omega_dot];

end

