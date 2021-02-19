clear;
clc;

%%
run CRJ330Config.m

%%
% Longitudinal

Long_data = CRJData.DimStableDer.Longitudinal;
save('long_dat.mat','-struct','Long_data');
load long_dat.mat

% Lateral

Lat_data = CRJData.DimStableDer.Lateral;
save('lat_dat.mat','-struct','Lat_data');
load lat_dat.mat

%%

I_xx = Iz/((Ix*Iz)-(Ixz*Ixz));
I_zz = Ix/((Ix*Iz)-(Ixz*Ixz));
I_xz = Ixz/((Ix*Iz)-(Ixz*Ixz));

%% Linearized longitudinal equation

A31 = (Mu + (Zu*Mw_dot)/(m-Zw_dot))/Iy;
A32 = (Mw + (Zw*Mw_dot)/(m-Zw_dot))/Iy;
A33 = (Mq + ((Zq+m*u0)*Mw_dot)/(m-Zw_dot))/Iy;
A34 = -(m*g*sin(theta0)*Mw_dot)/(Iy*(m-Zw_dot));
A_long = [Xu/m Xw/m Xq/m -g*cos(theta0); Zu/(m-Zw_dot) Zw/(m-Zw_dot) (Zq+m*u0)/(m-Zw_dot) -(m*g*sin(theta0))/(m-Zw_dot); A31 A32 A33 A34; 0 0 1 0];
B_long = [Xde/m Xdp/m; Zde/(m-Zw_dot) Zdp/(m-Zw_dot); (Mde+(Zde*Mw_dot)/(m-Zw_dot))/Iy (Mdp+(Zdp*Mw_dot)/(m-Zw_dot))/Iy; 0 0];

sys_long = ss(A_long, B_long, eye(4), zeros(4,2));
[eigvec_Along, eigval_Along] = eig(A_long);

%% Linearized lateral equation

A_lat = [Yv/m Yp/m (Yr/m)-u0 g*cos(theta0); I_xx*Lv+I_xz*Nv I_xx*Lp+I_xz*Np I_xx*Lr+I_xz*Nr 0; I_zz*Nv+I_xz*Lv I_zz*Np+I_xz*Lp I_zz*Nr+I_xz*Lr 0; 0 1 tan(theta0) 0];
B_lat = [Yda/m Ydr/m; I_xx*Lda+I_xz*Nda I_xx*Ldr+I_xz*Ndr; I_zz*Nda+I_xz*Lda I_zz*Ndr+I_xz*Ldr; 0 0];

sys_lat = ss(A_lat, B_lat, eye(4), zeros(4,2));
[eigvec_Alat, eigval_Alat] = eig(A_lat);

%% Dynamic Modes

% Longitudinal Modes
[eigvec_Along, eigval_Along] = eig(A_long);
% Eigen Vector
long_eigvec = eigvec_Along;
% Eigen Value
D_long = diag(eigval_Along);

% Lateral Modes
[eigvec_Alat, eigval_Alat] = eig(A_lat);
% Eigen Vector
lat_eigvec = eigvec_Alat;
% Eigen Value
D_lat = diag(eigval_Alat);

%% Short period Mode
eig_sp = D_long(1:2,1);
% The below commented code gives the same result
% A_sp = [A_long(2,2) A_long(2,3); A_long(3,2) A_long(3,3)];
% eig_sp = eig(A_sp);

%% Phugoid Mode
eig_ph = D_long(3:4,1);

%% Spiral Mode
eig_sprl = D_lat(4,1);

%% Rolling Mode
eig_rl = D_lat(1,1);

%% Dutch Roll Mode
eig_dr = D_lat(2:3,1);
% A_dutchroll = [A_lat(1,1) A_lat(1,3); A_lat(3,1) A_lat(3,3)];
% eig_dr = eig(A_dutchroll);

%% Altitude Hold
% The A and B matrix for the altitude hold system is developed below with
% the new state h (altitude)

A_alh = [[A_long; [sin(theta0) -cos(theta0) 0 u0*cos(theta0)]] [0;0;0;0;0]];
B_alh = [B_long; [0 0]];

%% Roll Control

A_r = [[[A_lat B_lat(:,1)]; [0 0 0 0 -20]; [0 0 0 -1 0]] [0;0;0;0;0;0]];
B_r = [0; 0; 0; 0; 20; 0];
G_r = [0; 0; 0; 0; 0; 1];
C_r = [0 1 0 0 0 0
    0 0 0 1 0 0
    0 0 0 0 0 1];
F_r = zeros(3,1);

sys_roll = ss(A_r,B_r,C_r,F_r);

Qr = 50*[100 0 0 0 0 0
    0 1 0 0 0 0
    0 0 1 0 0 0
    0 0 0 1 0 0
    0 0 0 0 1 0
    0 0 0 0 0 10];
Rr = 100;

[Kr,Sr,er] = lqr(sys_roll,Qr,Rr);
[K1r,S1r] = linsolve(C_r',Kr');
K_r = 10*K1r';

A_rf = A_r - B_r*K_r*C_r;
B_rf = G_r - B_r*K_r*F_r;

sys_clr = ss(A_rf, B_rf, C_r, zeros(3,1));
stepinfo('sys_clr');
figure;
stepplot(sys_clr);

t = 0:0.01:100;
[y,tt,xx1] = step(sys_clr,t);
z2r = [0 0 0 1 0 0]*xx1';

figure;
plot(t,z2r);
xlabel('time');
ylabel('phi');
grid on;

%% transfer function
 
[NumG1,Den] = ss2tf(A_r,B_r,C_r,zeros(3,1));
G11 = tf(NumG1(1,:),Den);
G21 = tf(NumG1(2,:),Den);
G31 = tf(NumG1(3,:),Den);
%% bode plot data: KG

figure;
bode(G11,'r',G21,'b',G31,'g');
grid on;
legend('G11','G21','G31');
htype = findobj(gcf,'type','line');
set(htype,'linewidth',2);

%% Singular Value

figure;
sigma(sys_roll);
grid on;

%%

for i = 1:1001
    w(i) = 0.1*i;
end

%% frequency response of KG

[mag1,phase1] = bode(G11,w);
[mag2,phase2] = bode(G21,w);
[mag3,phase3] = bode(G31,w);

mag_KG(1,:) = K_r(1,1)*mag1(1,1,:) + K_r(1,2)*mag2(1,1,:) + K_r(1,3)*mag3(1,1,:);

mag_KG_db1 = 20*log10(mag_KG(1,:));

%% bode plot data: Wind gust

L = 3.39;
sigma = 0.1;
a1 = 2*L*sigma^2;
a2 = 3*L^2;
a3 = L^2;
for i = 1:1001
    phi_w(i) = a1*(1+a2*w(i)^2)/(1+a3*(w(i)^2))^2;
end
phi_wdb = 20*log10(phi_w);

%% bode plot data: Model uncertainty

zeta = 0.3;
wn = 40;
M = tf([-1 -2*zeta*wn 0],[1 2*zeta*wn wn^2]);
[magM, phaseM] = bode(M,w);
magMdb(1,:) = -20*log10(magM(1,1,:));


%% frequency response and robustness simulation 

figure;
plot(log10(w), phi_wdb, log10(w), magMdb, log10(w), mag_KG_db1); 

legend('\phi_w(s)','1/M(s)', 'KG1');
xlabel('log10\omega');
title('Frequency Response');

%%%%%%%%%%%%-----------------------------%%%%%%%%%%%%%%%%%%%5