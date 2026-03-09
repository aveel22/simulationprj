function [u, state] = pid_controller(e, state, Kp, Ki, Kd, Ts)
    % pid_controller - Discrete PID for scalar or vector signals
    %
    % Inputs:
    %   e     - current error (scalar or vector)
    %   state - struct with fields 'integral' and 'prevError'
    %   Kp, Ki, Kd - PID gains (scalar or vector of same size as e)
    %   Ts    - sample time
    %
    % Outputs:
    %   u     - control signal (same size as e)
    %   state - updated state
    
    % Ensure state fields have correct shape
    if isempty(state.integral)
        state.integral = zeros(size(e));
    end
    if isempty(state.prevError)
        state.prevError = zeros(size(e));
    end
    
    % Update integral (elementwise)
    state.integral = state.integral + e .* Ts;
    
    % Derivative (elementwise)
    derivative = (e - state.prevError) ./ Ts;
    
    % PID law (elementwise)
    u = Kp .* e + Ki .* state.integral + Kd .* derivative;
    
    % Update memory
    state.prevError = e;
end


