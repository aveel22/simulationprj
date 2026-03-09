function [control_vector, current_error, integral_sum] = pid_controller(xh, xr, p, dt, integral_sum, prev_error)
    % PID Controller implementation for 3D vector control
    %
    % Inputs:
    %   xh - current value (3,1) or (1,3) vector
    %   xr - reference vector (3,1)
    %   p - parameter structure with p.ctrl.p, p.ctrl.i, p.ctrl.d as (3,1) vectors
    %   dt - time step for derivative and integral calculations
    %   integral_sum - accumulated integral error (3,1), [] for first call
    %   prev_error - previous error for derivative calculation (3,1), [] for first call
    %
    % Outputs:
    %   control_vector - PID control output (3,1)
    %   current_error - current error (3,1)
    %   integral_sum - updated integral sum for next iteration
    
    % Ensure proper shapes - convert to column vectors
    xh = reshape(xh, 3, 1);  % Convert to (3,1)
    xr = reshape(xr, 3, 1);  % Ensure (3,1) shape
    
    % Calculate current error
    current_error = xr - xh;
    
    % Proportional term
    proportional = p.ctrl.p .* current_error;
    
    % Integral term
    if isempty(integral_sum)
        integral_sum = zeros(3, 1);
    end
    integral_sum = integral_sum + current_error * dt;
    integral = p.ctrl.i .* integral_sum;
    
    % Derivative term
    if isempty(prev_error)
        derivative = zeros(3, 1);
    else
        derivative = p.ctrl.d .* (current_error - prev_error) / dt;
    end
    
    % Calculate control vector
    control_vector = proportional + integral + derivative;
end

% Example usage script (separate .m file or script section)
% 
% % Initialize parameters structure
% p.ctrl.p = [1.0; 1.2; 0.8];    % Proportional gains (3,1)
% p.ctrl.i = [0.1; 0.15; 0.05];  % Integral gains (3,1)
% p.ctrl.d = [0.05; 0.08; 0.03]; % Derivative gains (3,1)
% 
% % Example values
% xh = [1.0; 2.0; 1.5];          % Current position (3,1)
% xr = [2.0; 3.0; 2.0];          % Reference position (3,1)
% dt = 0.01;                     % Time step
% 
% % First call (no previous state)
% [control, error, integral_sum] = pid_controller(xh, xr, p, dt, [], []);
% 
% fprintf('Control vector:\n');
% disp(control);
% fprintf('Current error:\n');
% disp(error);
% fprintf('Integral sum:\n');
% disp(integral_sum);
% 
% % Subsequent calls would use the returned integral_sum and previous error
% % Example of next iteration:
% xh_new = [1.1; 2.1; 1.6];      % Updated current position
% [control_new, error_new, integral_sum_new] = pid_controller(...
%     xh_new, xr, p, dt, integral_sum, error);
% 
% fprintf('\n--- Next iteration ---\n');
% fprintf('Control vector:\n');
% disp(control_new);
% fprintf('Current error:\n');
% disp(error_new);

