function [A, B] = build_lsys(x, u, p)
%BUILD_LSYS Summary of this function goes here
%   Detailed explanation goes here
% x - state vector for linearization
% u - control vector
% p - parameters

%% Declare
A = zeros(size(x,1));
B = zeros(size(x,1), size(u,1));
r = x(1:3);                     % coordinates
v = x(4:6);                     % velocity
eul0 = x(7:9);                  % angles
w = x(10:12);                   % angular velocity
%% Unpack parameters
m = p.m;                % current mass
Mass = diag([m, m, m]);
Js = tensor_interpolation(m, p.mass, p.Inertia);
rcg = mass_center_interpolate(m, p.mass, p.rcgs);
[Hf, Hm] = ComponentForForces(p.geom, rcg);
Cv = p.Cv;
Cw = p.Cw;

A11 = zeros(3,3);
A12 = eye(3,3);
A21 = zeros(3,3);
A22 = -skew(w) + Mass \ (-Cv); % + Nrot * G; for linear system exclude gravity (to make system Linear)

At = [A11, A12;
      A21, A22];

A33 = zeros(3,3);
A43 = zeros(3,3);
A34 = Tw(eul0);
A44 = Js \ (-skew(w) * Js - Cw);
Ar = [A33, A34;
      A43, A44];

A = [At, zeros(6,6); 
    zeros(6,6), Ar];
if p.grav_on > 0
    A(4:6,7:9) = grav_lin(eul0, p.g);
end
%% Control
B11 = zeros(3, 9);
B21 = Mass \ Hf;
Bt = [B11;
      B21];
B31 = zeros(3,9);
B41 = Js \ Hm;
Br = [B31;
      B41];
B = [Bt; Br];
end

