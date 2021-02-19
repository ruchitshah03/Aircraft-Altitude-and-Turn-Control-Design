function [state_dot] = get_x_dot1(t, state, input)
%GET_X_DOT Summary of this function goes here
%   Detailed explanation goes here
% run CRJ330Config.m
% run data_load.m

load Data.mat

% D = create_D_struct();
% C = create_C_struct();
% Data = create_Data_struct();
% [C, Data] = set_crj_imperial(C, Data); 
% D = calc_D_struct(C, D, Data); 

u = state(1);
v = state(2);
w = state(3);

p = state(4);
q = state(5);
r = state(6);

phi = state(7);
theta = state(8);
psi = state(9); 

x = state(10);
y = state(11);
z = state(12); 

delta_e = input(1);
delta_a = input(2);
delta_r = input(3);
delta_p = input(4);

theta_e = theta0;
U_e = u0;

%%
X_stat = - m .* g .* sin(theta); 
Y_stat = m * g * sin(phi) * cos(theta);
Z_stat = m * g * cos(theta) * cos(phi);

w_dot = ((Z_stat - m .* g .* cos(theta_e) + Zu .* (u - U_e) + Zw .* w + Zq .* q + Zde .* delta_e)./m + q .* u - p .* v)./(1-Zw_dot ./ m);

X_aero = m .* g .* sin(theta_e) + Xu .* (u - U_e) + Xw .* w + Xde .* delta_e + Xdp .* delta_p;
Y_aero = Yv .* v + Ydr .* delta_r + Yp .* p + Yr .* r;
Z_aero = - m .* g .* cos(theta_e) + Zu .* (u - U_e) + Zw .* w + Zq .* q + Zw_dot .* w_dot + Zde .* delta_e;

X = X_stat + X_aero; 
Y = Y_stat + Y_aero;
Z = Z_stat + Z_aero;

L = Lv .* v + Lp .* p + Lr .* r + Lda .* delta_a + Ldr .* delta_r;
M = Mu .* (u - U_e) + Mw .* w + Mq .* q + Mw_dot .* w_dot + Mde .* delta_e;
N = Nv .* v + Np .* p + Nr .* r + Nda .* delta_a + Ndr .* delta_r;

u_dot = X./m + r .* v - w .* q;
v_dot = Y./m + w .* p - r .* u;

p_dot = (L + (N*Izx./Iz) + q.*r.*(Iy - Iz + ((Izx.^2)./Iz)) + p.*q.*(Ix-Iy+Iz).*(Izx./Iz)).*(Iz./(Ix.*Iz + (Izx.^2)));
q_dot = (M + Izx.*(p.^2 - r.^2) + (Iz-Ix) .* p .* r)./Iy;
% r_dot = (N + L .* Izx ./ Ix + p .* q .* (Ix - Iy - Izx.^2 ./ Ix) + q .* r .* (Iz-Iy-Ix) .* Izx ./ Ix) .* (Ix ./ (Ix .* Iz - Izx.^2));
r_dot = (N + Ixz.*(q.*r - p_dot) + (Ix-Iy).*p.*q)./Iz;

% Equation 1.119 page 34
phi_dot = p + sin(phi) .* tan(theta) .* q + cos(phi) .* tan(theta) .* r;
theta_dot = cos(phi) .* q - sin(phi) .* r;
psi_dot = sin(phi) .* sec(theta) .* q + cos(phi) .* sec(theta) .* r;    


% Equation 1.61 page 21
C_EB_1_1 = cos(theta) .* cos(psi); 
C_EB_1_2 = sin(phi) .* sin(theta) .* cos(psi) - cos(phi) .* sin(psi);
C_EB_1_3 = cos(phi) .* sin(theta) .* cos(psi) + sin(phi) .* sin(psi);

C_EB_2_1 = cos(phi) .* sin(psi);
C_EB_2_2 = sin(phi) .* sin(theta) .* sin(psi) + cos(phi) .* cos(psi);
C_EB_2_3 = cos(phi) .* sin(theta) .* sin(psi) - sin(phi) .* cos(psi);

C_EB_3_1 = -sin(theta); 
C_EB_3_2 = sin(phi) .* cos(theta);
C_EB_3_3 = cos(phi) .* cos(theta);

% Equation 1.66 page 22
x_dot = C_EB_1_1 .* u + C_EB_1_2 .* v + C_EB_1_3 .* w;
y_dot = C_EB_2_1 .* u + C_EB_2_2 .* v + C_EB_2_3 .* w;
z_dot = C_EB_3_1 .* u + C_EB_3_2 .* v + C_EB_3_3 .* w;


state_dot = [u_dot, v_dot, w_dot, ...
            p_dot, q_dot, r_dot, ...
            phi_dot, theta_dot, psi_dot, ...
            x_dot, y_dot, z_dot];
end

