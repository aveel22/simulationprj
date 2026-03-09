function xdot = MiniHopperODE_6dof(x, u, p)
    %MINIHOPPERODE_6DOF Summary of this function goes here
    %   Detailed explanation goes here

    v = x(4:6);  % velocity
    a = x(7:9);  % Euler's angles
    w = (10:12);  % angular velocity

    % Mass and Inertia matrices
    rcg = p.rc;
    Jo = p.J;
    m = p.m;
    % Damping matrices
    Cv = p.Cv;
    Cw = p.Cw;
    % Gravity
    G = p.gc .* N2B(a) * p.G;    % Get gravity in body frame
    Mass = eye(3) * m;

    % Control Forces
    F = u(1:3);                  % Forces
    M = u(4:6);                  % Moments

    % System of Equations
    xdot = zeros(x.size);
    xdot(1:3,1) = v;
    xdot(4:6,1) = -cross(w, v) + Mass \ (-Cv * v + F) + G;
    xdot(7:9,1) = E2W(a) * w;
    xdot(10:12,1) = Jo \ (-cross(w, Jo * w) - Cw * w + M);
end

