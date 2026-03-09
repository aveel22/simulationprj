function mrp = eul2MRP(eul, sequence)
    % Converts Euler angles to Modified Rodrigues Parameters (MRPs)
    % Input:
    %   eul - [roll, pitch, yaw] (ZYX default)
    %   sequence - Rotation sequence (e.g., 'ZYX', 'XYZ')
    % Output:
    %   mrp - 3x1 vector of MRPs

    % Convert Euler angles to quaternion
    q = eul2quat(eul, sequence);

    % Extract quaternion components
    q0 = q(1);
    qv = q(2:4);  % Vector part [q1; q2; q3]

    % Compute MRPs
    if q0 ~= 0
        mrp = qv / (1 + q0);
    else
        % If q0 = 0, use the shadow set to avoid singularity
        mrp = -qv / (1 - q0);
    end
end

