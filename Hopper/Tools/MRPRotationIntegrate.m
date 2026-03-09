function sigma_dot = MRPRotationIntegrate(sigma, omega)
%MRPROTATIONINTEGRATE Summary of this function goes here
%   Detailed explanation goes here
% Compute the norm of sigma
    sigma2 = sigma' * sigma;
    % Compute the G matrix
    G = (1 - sigma2) * eye(3) + 2 * skew(sigma) + 2 * sigma * sigma';
    % Compute sigma_dot (MRP derivative)
    sigma_dot = 0.25 * G * omega;
end