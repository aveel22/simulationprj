function Nrot = N2B(arg)
    % This is transformation matrix from Launch to Body Frame
    % Yaw(Z)->Pitch(Y)->Roll(X)
    % Unpack the input vector
    phi = arg(1);   % roll
    theta = arg(2); % pitch
    psi = arg(3);   % yaw
    Nrot = rotx(phi) * roty(theta) * rotz(psi);
end

