function [xdot, DCM, G] = NL_HopperDyn_v4_quat(x, u, Jo, Jdot, mdot, p, d)
%                                              x, uuu, Jo, Jdot, mdot, p, d
% NL_HOPPERDYN_1 Summary of this function goes here
%   Detailed explanation goes here
% x(14) - State Vector 6DoF: px,py,pz, vx,vy,vz, q0,q1,q2,q3, wx,wy,wz, m
% u(9)  - Control Vector: Fx,Fy,Fz,Mx,My,Mz
% p(20) - System Parameters: Jxx,Jyy,Jzz,Jxy,Jyz,Jxz,m,Cx,Cy,Cz,Cr,Cq,Cp,g,
%       - beta, gamma, rho, h, d, k_quat


%% State vector
% r = x(1:3);             % position
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
cx      = p(7); 
cy      = p(8); 
cz      = p(9); 
cr      = p(10);
cq      = p(11); 
cp      = p(12); 
gg      = p(13); 
% beta    = p(14); 
% gamma   = p(15); 
% rho     = p(16); 
k_quat  = p(20);      % for quaternion normalization
gc = p(21);           % turnON garvity

%% Mass and Inertia matrices
Mass = diag([x(14), x(14), x(14)]);
%% Damping matrices
Cv = diag([cx, cy, cz]);
Cw = diag([cr, cq, cp]);
% Cv = diag([0, 0, 0]);
% Cw = diag([0, 0, 0]);
%% DCM
DCM = quat2DCM(q);                              % inertial to body DCM
ang = quat2eul(q');                             % yaw pitch roll - > PYR
L2B = rotx(ang(3))*roty(ang(2))*rotz(ang(1)); % L2B - DCM = zeros(3)
%% Gravity
G = gc * DCM * [0; -gg; 0];                          % Get gravity in body frame    
%% System of Equations 
xdot = zeros(size(x,1), 1);                     % initialize
xdot(1:3, 1)    = v;
xdot(4:6, 1)    = -cross(w, v) + Mass \ (-Cv * v + F + Fd) + G;
xdot(7:10, 1)   = QuatRotationIntegrate(q, w, k_quat);
xdot(11:13, 1)  = Jo \ (-cross(w, Jo*w)  - Cw*w - Jdot * w + M + Md);
xdot(14, 1)        = -mdot;