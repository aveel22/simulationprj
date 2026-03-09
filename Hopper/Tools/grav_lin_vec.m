function D2B = grav_lin_vec(eul, g)
% Returns the derivative of N2B(eul) w.r.t. Euler angles,
    % so that D2B_from_euler(eul, g) = linearized gravity acceleration.
    %
    % Inputs:
    %   eul - [roll, pitch, yaw] Euler angles vector
    %   g   - gravity vector (3x1)
    %
    % Output:
    %   D2B - 3x3 matrix with columns [∂R/∂roll*g, ∂R/∂pitch*g, ∂R/∂yaw*g]
    
    roll = eul(1);
    pitch = eul(2);
    yaw = eul(3);
    
    Rx = rotx(roll);
    Ry = roty(pitch);
    Rz = rotz(yaw);
    
    dRx = drotx_da(roll);
    dRy = droty_da(pitch);
    dRz = drotz_da(yaw);
    
    % Derivative contributions
    dR_droll = dRx * Ry * Rz;
    dR_dpitch = Rx * dRy * Rz;
    dR_dyaw = Rx * Ry * dRz;
    
    % Apply to gravity vector
    c_roll = dR_droll * g;
    c_pitch = dR_dpitch * g;
    c_yaw = dR_dyaw * g;
    
    % Stack columns
    D2B = [c_roll, c_pitch, c_yaw];
end

