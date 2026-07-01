function [xdot, DCM_e, G] = NL_HopperDyn_v2(x, u, Jo, Jdot, m)
% NL_HOPPERDYN_1 Summary of this function goes here
%   Detailed explanation goes here
% x - State Vector 6DoF: px,py,pz, vx,vy,vz,roll,pitch,yaw,wx,wy,wz
% u - Control Vector: Fx,Fy,Fz,Mx,My,Mz
% Jo - tensor of inertia
% Jdot - dot tensor of inertia
% m - current mass


%% State vector
r = x(1:3);             % position
v = x(4:6);             % velocity
a = x(7:9);             % attitude angles
w = x(10:12);           % angular velocity
F = u(1:3);             % Forces
M = u(4:6);             % Moments
%% Damping matrices
% Cv = diag([cx, cy, cz]);
% Cw = diag([cr, cq, cp]);

%% Mass and Inertia matrices
% Jo
Mass = diag([m, m, m]);
if m < 29.99
    i = 10; % for debugging
end
% MassInv = inv(Mass);
% Joinv = inv(Jo);
%% DCM
DCM_e = Nw(a);
%% Gravity
G = Nw(a) * [-9.8065; 0; 0];                 % Get gravity in body frame        
%% System of Equations
xdot = zeros(size(x,1), 1);             % initialize
xdot(1:3, 1)    = v;
xdot(4:6, 1)    = -cross(w, v) + Mass \ (F);
xdot(7:9, 1)    = Tw(a) * w;
xdot(10:12, 1)  = Jo \ (-cross(w, Jo*w) - Jdot * w + M);