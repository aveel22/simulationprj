function sigma_dot = MRP_Integrate(sigma, omega)
%% Function for MRP kinematics with shadow switching
% it's better not use in symbolic
% Threshold for shadow MRP switch (commonly 1)
    sigma_threshold = 1;
    
    % Compute the squared norm of sigma
    sigma2 = sigma' * sigma;

    % If norm exceeds threshold, switch to shadow MRP
    if sigma2 > sigma_threshold
        sigma = -sigma / sigma2;  % Shadow MRP transformation
        sigma2 = sigma' * sigma;  % Update norm after switching
    end
   
    % Compute the G matrix
    G = (1 - sigma2) * eye(3) + 2*skew(sigma) + 2 * sigma * sigma';
    
    % Compute sigma_dot (MRP derivative)
    sigma_dot = 0.25 * G * omega;
end

