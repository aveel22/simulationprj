function euler_continuous = quat_to_euler_continuous(q, eul_prev)
    % Converts quaternion to continuous Euler angles avoiding -180° to 180° jumps.
    % Inputs:
    %   q - Quaternion [q0; q1; q2; q3] (4x1 vector)
    %   prev_euler - Previous Euler angles [phi; theta; psi] (3x1 vector)
    % Outputs:
    %   euler_continuous - Continuous Euler angles [phi; theta; psi] (3x1 vector)
    
    % Convert quaternion to Euler angles (assuming ZYX convention)
    euler = quat2eul(q', 'ZYX')';  % Convert and transpose to column vector
    
    % Handle angle discontinuities by unwrapping
    euler_continuous = unwrap_angles(euler, eul_prev);
end

