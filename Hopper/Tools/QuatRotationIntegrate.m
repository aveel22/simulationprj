function qdot = QuatRotationIntegrate(quat, w, k_quat)
    %% Rotation Integration Using Quaternions
    % quat = [q0; q1; q2; q3] - quaternion
    % w = [wx; wy; wz] - angular velocity in body frame
    % k_quat - correction factor for unit norm maintenance (default: 1)

    % Default value for k_quat if not provided
    if nargin < 3
        k_quat = 1;
    end

    % Quaternion Multiplication Matrix (for integrating rotation)
    QM = [ -quat(2), -quat(3), -quat(4);
            quat(1), -quat(4),  quat(3);
            quat(4),  quat(1), -quat(2);
           -quat(3),  quat(2),  quat(1)];
    
    % Compute the norm correction term
    qnorm = (1 - k_quat * (quat' * quat));  % Ensuring unit norm correction

    % Compute the quaternion derivative
    qdot = 0.5 * QM * w + qnorm * quat;