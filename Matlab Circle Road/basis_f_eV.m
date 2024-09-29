function basis_V = basis_f_eV(X) 
%% Get the state
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

%% 
basis = [V_y*omega, V_x^2 - V_L^2, V_x*omega - V_L*omega_L, V_y/V_x, omega_L/V_L - omega/V_x]';
end