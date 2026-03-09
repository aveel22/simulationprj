function res = Tw(arg)
    % Unpack the input vector
    % Rocket transformatio
    % Pitch(Z) -> Yaw(Y) -> Roll(Z)
    % Use of matrix
    % [phi; theta; psi] = TW(arg) * [wx; wy; wz]
    phi = arg(1);   % roll
    theta = arg(2); % pitch
    psi = arg(3);   % yaw
    
    % Check for singularity (theta = pi/2 or -pi/2 leads to division by zero)
    if cos(psi) == 0
        error('Theta must not be +/- pi/2 to avoid division by zero.');
    end
    
    % Compute the transformation matrix
    res = (1 / cos(psi)) * [
        cos(psi), sin(psi) * sin(phi), sin(psi) * cos(phi);
        0, sin(phi), cos(phi); 
        0, cos(psi) * cos(phi), -cos(psi) * sin(phi)
    ];
end