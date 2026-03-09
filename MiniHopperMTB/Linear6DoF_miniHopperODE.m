function [A, B] = Linear6DoF_miniHopperODE(x, u, p)
    %LINEAR6DOF_MINIHOPPERODE Summary of this function goes here
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
    
    if p.inv.roll
        eul0(1) = -eul0(1);
    end
    if p.inv.pitch
        eul0(2) = -eul0(2);
    end
    if p.inv.yaw
        eul0(3) = -eul0(3);
    end
    %% Unpack parameters
    m = p.m;                % current mass
    Mass = diag([m, m, m]);
    Js = p.J;
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
    A34 = E2W(eul0);
    A44 = Js \ (-skew(w) * Js - Cw);
    Ar = [A33, A34;
          A43, A44];
    
    A = [At, zeros(6,6); 
        zeros(6,6), Ar];
    %     if p.grav_on > 0
    %         A(4:6,7:9) = grav_lin(eul0(1), eul0(2), eul0(3), p.g);
    %     end
    %% Control
    B11 = zeros(3, 3);
    B21 = Mass \ eye(3);
    Bt = [B11;
          B21];
    B31 = zeros(3,3);
    B41 = Js \ eye(3);
    Br = [B31;
          B41];
    B = [Bt, zeros(6,3); zeros(6,3), Br];
end
