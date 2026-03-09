function dcm = mrp2DCM(mrp)
%MRP2DCM Summary of this function goes here
%   Detailed explanation goes here
% Compute DCM from MRPs
    mrp2 = mrp' * mrp;
    I = eye(3);
    
    % DCM using MRPs
    dcm = I + (8 / (1 + mrp2)) * (skew(mrp) * skew(mrp) - 0.5 * skew(mrp));
end

