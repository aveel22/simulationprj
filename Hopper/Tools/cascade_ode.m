function dx = cascade_ode(t, xs, u, T, n)
    % cascade_ODE - Cascade of 1st-order filters
    %
    % Inputs:
    %   t  - time (not used, but needed for ODE solvers)
    %   xs - state vector (n×1)
    %   u  - input signal (scalar)
    %   T  - time constant
    %   n  - number of cascaded stages
    %
    % Output:
    %   dx - derivative of states (n×1)

    dx = zeros(size(xs));  % same shape as xs

    % Adjust time constant if last state > input
    if xs(end) > u
        T = T / 4;
    end

    % First stage
    dx(1) = (u - xs(1)) / T;

    % Remaining cascaded stages
    for i = 2:n
        dx(i) = (xs(i-1) - xs(i)) / T;
    end
end


