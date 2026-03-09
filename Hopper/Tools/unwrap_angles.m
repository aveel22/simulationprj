function euler_smooth = unwrap_angles(euler, eul_prev)
    % Ensures Euler angles remain continuous by unwrapping jumps.
    % Inputs:
    %   euler - Current Euler angles [phi; theta; psi] (3x1 vector)
    %   prev_euler - Previous Euler angles [phi; theta; psi] (3x1 vector)
    % Outputs:
    %   euler_smooth - Smoothed Euler angles (3x1 vector)

    euler_smooth = euler;
    
    % Check each angle and correct discontinuities
    for i = 1:3
        while euler_smooth(i) - eul_prev(i) > 180
            euler_smooth(i) = euler_smooth(i) - 360;
        end
        while euler_smooth(i) - eul_prev(i) < -180
            euler_smooth(i) = euler_smooth(i) + 360;
        end
    end
end
