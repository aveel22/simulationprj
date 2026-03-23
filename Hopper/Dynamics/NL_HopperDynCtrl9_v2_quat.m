function [xdot, G, DCM] = NL_HopperDynCtrl9_v2_quat(x, u, Jo, Jdot, m, p)
%% Used for 6DoF_v3_quat/LQR_6DoF_v3.slx
% NL_HOPPERDYN_1 Summary of this function goes here
%   Detailed explanation goes here
% x(14) - State Vector 6DoF: px,py,pz, vx,vy,vz, q0,q1,q2,q3, wx,wy,wz, m
% u(9)  - Control Vector: F1,F2,F3,F4,F5,F6,F7,F8,F9
% p(20) - System Parameters: Jxx,Jyy,Jzz,Jxy,Jyz,Jxz,m,Cx,Cy,Cz,Cr,Cq,Cp,g,
%       - beta, gamma, rho, h, d, k_quat
%% State vector
r = x(1:3);             % position
v = x(4:6);             % velocity
q = x(7:10);            % quaternion
w = x(11:13);           % angular velocity


%% Parameters
cx      = p(1); 
cy      = p(2); 
cz      = p(3); 
cr      = p(4);
cq      = p(5); 
cp      = p(6); 
gg      = p(7); 
% rho     = p(8); 
h       = p(9);    % roll levers
f       = p(10);    % yaw-pitch levers
k_quat  = p(11);
%% Mass parameter
Mass = diag([m, m, m]);

%% Levers and forces
sq2 = 1 / sqrt(2);
% Force distrubution matrix
%     F1  F2     F3    F4    F5 F6 F7 F8 F9
Hf = [1,   0,     0,    0,    0, 1, 1, 1, 1;    % Fx
      0, sq2,   sq2, -sq2, -sq2, 0, 0, 0, 0;    % Fy
      0, -sq2,  sq2,  sq2, -sq2, 0, 0, 0, 0];   % Fz

% Lever matrix
%     F1 F2 F3  F4 F5 F6 F7  F8 F9
Hm = [0, -h, h, -h, h, 0, 0,  0, 0;
      0,  0, 0,  0, 0, f, 0, -f, 0;
      0,  0, 0,  0, 0, 0, -f, 0, f];

F = Hf * u;             % Forces
M = Hm * u;             % Moments
%% Damping matrices
Cv = diag([cx, cy, cz]);
Cw = diag([cr, cq, cp]);

%% DCM
DCM = quat2DCM(q);                              % inertial to body DCM
ang = quat2eul(q');                             % yaw pitch roll
L2B = rotx(ang(3))*roty(ang(2))*rotz(ang(1));   % L2B
%% Gravity
G = DCM * [0; 0; gg];                           % Get gravity in body frame    
%% System of Equations
xdot = zeros(size(x,1), 1);             % initialize
xdot(1:3, 1)    = v;
xdot(4:6, 1)    = -cross(w, v) + Mass \ (-Cv * v + F); % + G;
xdot(7:10, 1)   = QuatRotationIntegrate(q, w, k_quat);
xdot(11:13, 1)  = Jo \ (-cross(w, Jo*w)  - Cw*w - Jdot * w + M);