function eul = MRP2eul(mrp, sequence)
    % Converts Modified Rodrigues Parameters (MRPs) to Euler angles
    % Input:
    %   mrp - 3x1 MRP vector
    %   sequence - Rotation sequence (e.g., 'ZYX')
    % Output:
    %   eul - [roll, pitch, yaw] in radians

    % Compute equivalent quaternion
    q0 = (1 - norm(mrp)^2) / (1 + norm(mrp)^2);
    qv = (2 * mrp) / (1 + norm(mrp)^2);
    q = [q0; qv];

    % Convert to Euler angles
    eul = quat2eul(q', sequence);
end

