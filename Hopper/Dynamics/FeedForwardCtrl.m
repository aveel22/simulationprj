function [uuu] = FeedForwardCtrl(x, u, param, xddot)
    %FEEDFORWARDCTRL Summary of this function goes here
    %   Detailed explanation goes here
    %% Control
    uuu = zeros(size(u));
    if nargmin < 4
        xddot = zeros(6,1);
    end
    %% State vector
    % r = x(1:3);           % position
    v = x(4:6);             % velocity
    a = x(7:9);             % Euler's angles
    w = x(10:12);           % angular velocity
    m = x(13);              % current mass
    
    %% Mass and Inertia matrices    
    %     rcg = mass_center_interpolate(x(13), param.mass, param.rcgs);
    %     [Hf, Hm] = ComponentForForces(param.geom, rcg);
    %% Control Forces
    %     uuu = [Hf; Hm] * u; 

    Jo = tensor_interpolation(m, param.mass, param.Inertia);
    Jdot = param.Jdot_mlt * param.mdot;               % debug
    
    Mass = diag([m, m, m]);
    %% Damping matrices
    Cv = param.Cv;
    Cw = param.Cw;
    %% DCM
    DCM = Nw(a);                            % inertial to body DCM
    %% Gravity
    G = param.gc * DCM * [0; -param.g; 0];  % Get gravity in body frame
    %% FeedForward Control based on dynamic
    F = xddot(1:3) + cross(w, v) + Mass \ (Cv * v) - G;
    M  = xddot(4:6) + Jo \ (cross(w, Jo*w)  + Cw*w + Jdot * w);
    uuu = [F; M];
end

