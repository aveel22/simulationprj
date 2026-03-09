function Lrot = Lw(arg)
    % This is transformation matrix from Launch to Body Frame
    % Pitch(Z)->Yaw(Y)->Roll(X)
    % Unpack the input vector
    % Linearized DCM    
    % Appazov, p. 45, (2.27)
    phi = arg(1);   % roll
    theta = arg(2); % pitch
    psi = arg(3);   % yaw

    Lrot = [
            cos(theta),     sin(theta),                     -psi;
           -sin(theta),     cos(theta),                      phi;
           phi*sin(theta), -phi*cos(theta)+psi*sin(theta),   1
    ];
end