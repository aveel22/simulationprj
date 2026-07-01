function xdot = NL_HopperDyn_1(x, u, p, d)
% NL_HOPPERDYN_1 Summary of this function goes here
%   Detailed explanation goes here
% x - State Vector 6DoF: px,py,pz, vx,vy,vz,roll,pitch,yaw,wx,wy,wz
% u - Control Vector: Fx,Fy,Fz,Mx,My,Mz
% p - System Parameters: Jxx,Jyy,Jzz,Jxy,Jyz,Jxz,m,Cx,Cy,Cz,Cr,Cq,Cp,g
% d - Disturbances: dFx, dFy, dFz, dMx, dMy, dMz


r = x(1:3);     % position
v = x(4:6);     % velocity
a = x(7:9);     % attitude angles
w = x(10:12);   % angular velocity
F = u(1:3);    % Forces
M = u(4:6);    % Moments

% Disturbances
Fd = d(1:3);        % Disturbance Forces
Md = d(4:6);        % Disturbance Moments

% Unpack parameters
Jxx  = p(1);
Jyy  = p(2);
Jzz  = p(3);
Jxy  = p(4);
Jyz  = p(5);
Jxz  = p(6);
m    = p(7);
cx   = p(8);
cy   = p(9);
cz   = p(10);
cr   = p(11);
cq   = p(12);
cp   = p(13);
g    = p(14);
beta  = p(15);
gamma = p(16);
ro    = 1.2255;     % density: we can put it into parameters later
b = beta * F(1);    % mdot

% Damping matrices
Cv = diag([cx, cy, cz]);
Cw = diag([cr, cq, cp]);

% Inertia matrices
Jo = [Jxx, -Jxy, -Jxz;
     -Jxy, Jyy, -Jyz;
     -Jxz, -Jyz, Jzz];
Joinv = inv(Jo);
Mass = diag([m, m, m]);
Minv = inv(Mass);

% Gravity
G = Nw(a) * [0; -g; 0]; % Get gravity in body frame

% step
xdot = zeros(size(x,1), 1); % initialize
% System of Equations
xdot(1:3, 1)    = v;
xdot(4:6, 1)    = -cross(w, v) + Minv * (-Cv * v + F) + Fd;
xdot(7:9, 1)    = Tw(a) * w;
xdot(10:12, 1)  = Joinv * (-cross(w, Jo*w) - Cw*w + M) + Md;
