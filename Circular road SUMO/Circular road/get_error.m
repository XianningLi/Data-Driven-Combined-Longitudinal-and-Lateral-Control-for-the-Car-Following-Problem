function e = get_error(X)

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
radius =  51.6;
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
e1 = x_L + s*sin(phi_L) - x - d*cos(phi);
e2 = y_L - s*cos(phi_L) - y - d*sin(phi);
e3 = phi_L - (phi + gamma);
e4 = V_L - V_x;
e5 = -V_y;
e6 = omega_L - omega;
z1 =  cos(phi+gamma)*e1 + sin(phi+gamma)*e2;
z2 = -sin(phi+gamma)*e1 + cos(phi+gamma)*e2;

e = [z1;z2;e3;e4;e5;e6];

end