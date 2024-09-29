function basis_ue = get_basis_ue(X) 
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

e = get_error(X);

%%
basis_ue = [];
Vec1 = e';

CS0 = 1;
CS1 = [(1-cos(e(3))); sin(e(3))];
CS2 = [cos(e(3))^2; cos(e(3))*sin(e(3))];
CS3 = [cos(e(3))^3; cos(e(3))^2*sin(e(3))];
CS4 = [cos(e(3))^4; cos(e(3))^3*sin(e(3))];
CS5 = [cos(e(3))^5; cos(e(3))^4*sin(e(3))];
CS6 = [cos(e(3))^6; cos(e(3))^5*sin(e(3))];

% combs = nmultichoosek(Vec1, 5);
% Poly5 = combs(:,1).*combs(:,2).*combs(:,3).* combs(:,4).*combs(:,5);
% basis_ue = [basis_ue;Poly5];
% 
% combs = nmultichoosek(Vec1, 4);
% Poly4 = combs(:,1).*combs(:,2).*combs(:,3).* combs(:,4);
% basis_ue = [basis_ue;Poly4];

% combs = nmultichoosek(Vec1, 3);
% Poly3 = combs(:,1).*combs(:,2).*combs(:,3);
% basis_ue = [basis_ue; Poly3];

combs = nmultichoosek(Vec1, 2);
Poly2 = combs(:,1).*combs(:,2);
basis_ue = [basis_ue; Poly2];

combs = nmultichoosek(Vec1, 1);
Poly1 = combs(:,1);
basis_ue = [basis_ue; Poly1; kron(Poly1,CS1);];

% basis_ue = [basis_ue; CS1; CS2; CS3; 1]; 
basis_ue = [basis_ue; e(4)/(V_x); e(5)/(V_x); e(6)/(V_x)];


% basis_ue = [basis_ue; basis_fev'];
%%
% e123 = [e1, e2, e3];
% e456 = [e4, e(5), e6];
% 
% 
% basis_fep = [V_L*cos(phi_L), cos(phi_L)*omega_L, V_L*cos(phi_L-e3), V_L*sin(phi_L-e3), cos(phi_L-e3), sin(phi_L-e3),V_L*sin(phi_L), sin(phi_L)*omega_L];
% basis_fepdot = [V_L*sin(phi_L)*omega_L, sin(phi_L)*omega_L^2, V_L*sin(phi_L-e3)*(omega_L-e6), V_L*cos(phi_L-e3)*(omega_L-e6), sin(phi_L-e3)*(omega_L-e6), cos(phi_L-e3)*(omega_L-e6),V_L*cos(phi_L)*omega_L, cos(phi_L)*omega_L^2];
% basis_gep_fep = [basis_fep*sin(phi_L-e3),basis_fep*cos(phi_L-e3),basis_fep*sin(2*(phi_L-e3)),basis_fep*cos(2*(phi_L-e3))];
% basis_gep_e123 = [e123, e123*cos(phi_L-e3), e123*sin(phi_L-e3), e123*cos(2*(phi_L-e3)), e123*sin(2*(phi_L-e3))];
% basis_gep_e456 = [e456, e456*cos(phi_L-e3), e456*sin(phi_L-e3), e456*cos(2*(phi_L-e3)), e456*sin(2*(phi_L-e3))];
% basis_gepdot_fep = basis_gep_fep*(omega_L-e6);
% basis_gepdot_e123 = basis_gep_e123*(omega_L-e6);
% basis_gep_fepdot = [basis_fepdot*sin(phi_L-e3), basis_fepdot*cos(phi_L-e3),basis_fepdot*sin(2*(phi_L-e3)), basis_fepdot*cos(2*(phi_L-e3))];
% basis_gep_e123dot = e456; 
% basis_fev = [e(5)*omega_L, e(5)*e6, e4^2, e4*V_L, e4*e6, e4*omega_L, e6*V_L, e(5)/(V_L-e4), omega_L/V_L, (omega_L - e6)/(V_L - e4), e(5)/(V_L - e4)];
% 
% basis_ue = [basis_fep'; basis_fepdot'; basis_gep_e456'; basis_gep_fep'; basis_gep_e123'; basis_gepdot_fep'; basis_gepdot_e123'; basis_gep_fepdot'; basis_gep_e123dot'; basis_fev'];
end