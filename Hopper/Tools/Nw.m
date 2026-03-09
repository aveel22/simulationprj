function Nrot = Nw(arg)
    % This is transformation matrix from Launch to Body Frame
    % Pitch(Z)->Yaw(Y)->Roll(X)
    % Unpack the input vector
    phi = arg(1);   % roll
    theta = arg(2); % pitch
    psi = arg(3);   % yaw

    % rotation around X axis
    rotx = [1, 0, 0;
            0, cos(phi), sin(phi);
            0, -sin(phi), cos(phi)];
    
    % rotation around Y axis
    roty = [cos(psi), 0, -sin(psi);
            0,        1,    0;
            sin(psi), 0,    cos(psi)];
    
    % rotation around Z axis
    rotz = [cos(theta), sin(theta), 0;
            -sin(theta), cos(theta), 0;
              0,          0,         1];
    Nrot = rotx * roty * rotz;
end