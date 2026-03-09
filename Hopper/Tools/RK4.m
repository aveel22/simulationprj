function x_next = RK4(func, t, x0, u, param)
% RK4 Integrates a nonlinear ODE using 4th-order Runge-Kutta method
%
%   x_next = RK4(func, t, x0, u, param)
%
% Inputs:
%   func  - function handle for ODE: dxdt = func(x, u, param)
%   t     - [t0 t1], initial and final time (e.g., [0 Ts])
%   x0    - initial state vector
%   u     - control input
%   param - system parameters
%
% Output:
%   x_next - next state after integration over [t0, t1]

    dt = t(2) - t(1);  % time step

    k1 = func(x0, u, param);
    k2 = func(x0 + 0.5 * dt * k1, u, param);
    k3 = func(x0 + 0.5 * dt * k2, u, param);
    k4 = func(x0 + dt * k3, u, param);

    x_next = x0 + (dt / 6) * (k1 + 2*k2 + 2*k3 + k4);
end


