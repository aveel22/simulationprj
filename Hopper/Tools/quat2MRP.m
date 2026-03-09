function mrp = quat2MRP(q)
%QUAT2MRP Summary of this function goes here
%   Detailed explanation goes here
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

