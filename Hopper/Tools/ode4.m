function x = ode4(ode, tspan, x0, u, p)
    % ODE4 Summary of this function goes here
    %   Detailed explanation goes here
    x = zeros(size(x0,1), size(tspan, 2));
    x(:,1) = x0;
    N = length(tspan)-1;
    for i = 1:N
        x(:, i+1) = RK4(ode, x(:, i), u(:, i), p);
    end
end

