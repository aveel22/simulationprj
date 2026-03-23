function [xdot, DCM, G] = NL_HopperDyn_v3_quat(x, u, Jo, Jdot, m, p, d)
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
%% Control Forces
F = u(1:3);             % Forces
M = u(4:6);             % Moments
%% Disturbances
Fd = d(1:3);            % Disturbance Forces
Md = d(4:6);            % Disturbance Moments
%% Parameters
% Jxx     = p(1); 
% Jyy     = p(2); 
% Jzz     = p(3); 
% Jxy     = p(4); 
% Jyz     = p(5); 
% Jxz     = p(6); 
% m       = p(7); 
cx      = p(8); 
cy      = p(8); 
cz      = p(8); 
cr      = p(11);
cq      = p(11); 
cp      = p(11); 
% gg      = p(14); 
% beta    = p(15); 
% % gamma   = p(16); 
% % rho     = p(17); 
% h       = p(18);    % roll levers
% f       = p(19);    % yaw-pitch levers
% k_quat  = p(20);      % for quaternion normalization
k_quat = 1;

%% Mass and Inertia matrices
Mass = diag([m, m, m]);
%% Levers and forces
sq2 = 1 / sqrt(2);
% Force distrubution matrix
%     F1  F2     F3    F4    F5 F6 F7 F8 F9
% Hf = [1,   0,     0,    0,    0, 1, 1, 1, 1;    % Fx
%       0, sq2,   sq2, -sq2, -sq2, 0, 0, 0, 0;    % Fy
%       0, -sq2,  sq2,  sq2, -sq2, 0, 0, 0, 0];   % Fz
% 
% % Lever matrix
% %     F1 F2 F3  F4 F5 F6 F7  F8 F9
% Hm = [0, -h, h, -h, h, 0, 0,  0, 0;
%       0,  0, 0,  0, 0, f, 0, -f, 0;
%       0,  0, 0,  0, 0, 0, -f, 0, f];
%% Damping matrices
Cv = diag([cx, cy, cz]);
Cw = diag([cr, cq, cp]);
% Cv = diag([0, 0, 0]);
% Cw = diag([0, 0, 0]);
%% DCM
DCM = quat2DCM(q);                              % inertial to body DCM
ang = quat2eul(q');                             % yaw pitch roll
% L2B = rotx(ang(3))*roty(ang(2))*rotz(ang(1));   % L2B
%% Gravity
G = DCM * [0; -9.8065; 0];                      % Get gravity in body frame    
% b = beta * u(1);                              % mdot
%% System of Equations 
xdot = zeros(size(x,1), 1);             % initialize
xdot(1:3, 1)    = v;
xdot(4:6, 1)    = -cross(w, v) + Mass \ (-Cv * v + F + Fd) + G;
xdot(7:10, 1)   = QuatRotationIntegrate(q, w, k_quat);
xdot(11:13, 1)  = Jo \ (-cross(w, Jo*w)  - Cw*w - Jdot * w + M + Md);