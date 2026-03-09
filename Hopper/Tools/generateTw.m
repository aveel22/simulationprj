function Tw = generateTw(arg)
    % Unpack the input vector
    psi = arg(1);
    theta = arg(2);
    phi = arg(3);
    
    % Check for singularity (theta = pi/2 or -pi/2 leads to division by zero)
    if cos(theta) == 0
        error('Theta must not be +/- pi/2 to avoid division by zero.');
    end
    
    % Compute the transformation matrix
    Tw = (1 / cos(theta)) * [
        0, sin(phi), cos(phi);
        0, cos(theta) * cos(phi), -cos(theta) * sin(phi);
        cos(theta), sin(theta) * sin(phi), sin(theta) * cos(phi)
    ];
end