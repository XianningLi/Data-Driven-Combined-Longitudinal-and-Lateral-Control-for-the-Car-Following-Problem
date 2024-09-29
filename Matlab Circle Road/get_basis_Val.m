function basis_Val = get_basis_Val(X) 
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
basis_Val = [];
Vec1 = e';

CS0 = 1;
CS1 = [(1-cos(e(3))); sin(e(3))];
CS2 = [(1-cos(e(3)))^2; cos(e(3))*sin(e(3))];
CS3 = [cos(e(3))^3; cos(e(3))^2*sin(e(3))];
CS4 = [cos(e(3))^4; cos(e(3))^3*sin(e(3))];
CS5 = [cos(e(3))^5; cos(e(3))^4*sin(e(3))];
CS6 = [cos(e(3))^6; cos(e(3))^5*sin(e(3))];

% combs = nmultichoosek(Vec1, 6);
% Poly6 = combs(:,1).*combs(:,2).*combs(:,3).* combs(:,4).*combs(:,5).* combs(:,6);
% basis_Val = [basis_Val;Poly6];
% 
% combs = nmultichoosek(Vec1, 5);
% Poly5 = combs(:,1).*combs(:,2).*combs(:,3).* combs(:,4).*combs(:,5);
% basis_Val = [basis_Val;Poly5];

% combs = nmultichoosek(Vec1, 4);
% Poly4 = combs(:,1).*combs(:,2).*combs(:,3).* combs(:,4);
% basis_Val = [basis_Val;Poly4];

combs = nmultichoosek(Vec1, 3);
Poly3 = combs(:,1).*combs(:,2).*combs(:,3);
basis_Val = [basis_Val; Poly3];

combs = nmultichoosek(Vec1, 2);
Poly2 = combs(:,1).*combs(:,2);
basis_Val = [basis_Val; Poly2; kron(Poly2,CS1)];

combs = nmultichoosek(Vec1, 1);
Poly1 = combs(:,1);
basis_Val = [basis_Val; Poly1; kron(Poly1,CS1)]; 

% basis_Val = [basis_Val; CS1; CS2; CS3]; 

basis_Val = [basis_Val; e(4)^2/(V_L-e(4)); e(5)^2/(V_L-e(4));e(6)^2/(V_L-e(4)); e(4)*e(5)/(V_L-e(4)); e(4)*e(6)/(V_L-e(4));e(5)*e(6)/(V_L-e(4))]; 

% combs = nmultichoosek(Vec1, 1);
% Poly1 = combs(:,1);
% basis_Val = [basis_Val; Poly1; kron(Poly1,CS1); kron(Poly1,CS2); kron(Poly1,CS3)]; 
% 
% basis_Val = [basis_Val; CS1; CS2; CS3]; 

%%
% e123 = [e1, e2, e3];
% e456 = [e4, e5, e6];
% 
% basis_fep = [V_L*cos(phi_L), cos(phi_L)*omega_L, V_L*cos(phi_L-e3), V_L*sin(phi_L-e3), cos(phi_L-e3), sin(phi_L-e3),V_L*sin(phi_L), sin(phi_L)*omega_L];
% basis_fepdot = [V_L*sin(phi_L)*omega_L, sin(phi_L)*omega_L^2, V_L*sin(phi_L-e3)*(omega_L-e6), V_L*cos(phi_L-e3)*(omega_L-e6), sin(phi_L-e3)*(omega_L-e6), cos(phi_L-e3)*(omega_L-e6),V_L*cos(phi_L)*omega_L, cos(phi_L)*omega_L^2];
% basis_gep_fep = [basis_fep*sin(phi_L-e3),ba sis_fep*cos(phi_L-e3),basis_fep*sin(2*(phi_L-e3)),basis_fep*cos(2*(phi_L-e3))];
% basis_gep_e123 = [e123, e123*cos(phi_L-e3),e123*sin(phi_L-e3), e123*cos(2*(phi_L-e3)), e123*sin(2*(phi_L-e3))];
% basis_gep_e456 = [e456, e456*cos(phi_L-e3), e456*sin(phi_L-e3), e456*cos(2*(phi_L-e3)), e456*sin(2*(phi_L-e3))];
% basis_gepdot_fep = basis_gep_fep*(omega_L-e6);
% basis_gepdot_e123 = basis_gep_e123*(omega_L-e6);
% basis_gep_fepdot = [basis_fepdot*sin(phi_L-e3), basis_fepdot*cos(phi_L-e3),basis_fepdot*sin(2*(phi_L-e3)), basis_fepdot*cos(2*(phi_L-e3))];
% basis_gep_e123dot = e456; 
% basis_fev = [e5*omega_L, e5*e6, e4^2, e4*V_L, e4*e6, e4*omega_L, e6*V_L, e5/(V_L-e4), omega_L/V_L, (omega_L - e6)/(V_L - e4), e5/(V_L - e4)];
% 
% 
% basis_ue = [basis_fep, basis_fepdot, basis_gep_e456, basis_gep_fep, basis_gep_e123, basis_gepdot_fep, basis_gepdot_e123, basis_gep_fepdot, basis_gep_e123dot, basis_fev];
% 
% basis_e1 = basis_ue*e1;
% basis_e2 = basis_ue*e2;
% basis_e3 = basis_ue*e3;
% basis_e4 = basis_ue*e4;
% basis_e5 = basis_ue*e5;
% basis_e6 = basis_ue*e6;
% 
% basis_Val = [basis_e1, basis_e2, basis_e3, basis_e4, basis_e5, basis_e6]';
% Vec = [e1, e2, e3, e4, e5, e6];
% combs = nmultichoosek(Vec, 4);
% basis_Val = [basis_Val; combs(:,1).*combs(:,2).*combs(:,3).* combs(:,4)];
% combs = nmultichoosek(Vec, 3);
% basis_Val = [basis_Val; combs(:,1).*combs(:,2).*combs(:,3)];
end