function [xhat_new, P_new] = ekf_d(xhat, P, u, y, ekf, p)
    % EKF_UPDATE_D Discrete Extended Kalman Filter (EKF) update step.
    %
    % INPUTS:
    %   xhat    - Current state estimate (nx x 1)
    %   P       - Current error covariance (nx x nx)
    %   u       - Control input (nu x 1)
    %   y       - Measurement vector (ny x 1)
    %   Ad      - Discrete state transition matrix (nx x nx)
    %   Bd      - Discrete control input matrix (nx x nu)
    %   C_meas  - Measurement matrix (ny x nx)
    %   Q_proc  - Process noise covariance (nx x nx)
    %   R_meas  - Measurement noise covariance (ny x ny)
    %
    % OUTPUTS:
    %   xhat_new - Updated state estimate (nx x 1)
    %   P_new    - Updated error covariance (nx x nx)
    
    C_meas = ekf.Cd; Q_proc = ekf.Qd; R_meas = ekf.Rd;
    Ad = p.plant.Ad; Bd = p.plant.Bd;

    % --- Prediction Step ---
    x_pred = Ad * xhat + Bd * u;
    P_pred = Ad * P * Ad' + Q_proc;
    
    % --- Update Step ---
    S = C_meas * P_pred * C_meas' + R_meas;
    K_kalman = P_pred * C_meas' / S;                    % Equivalent to P_pred * C_meas' * inv(S)
    y_residual = y - C_meas * x_pred;
    
    xhat_new = x_pred + K_kalman * y_residual;
    P_new = (eye(length(xhat)) - K_kalman * C_meas) * P_pred;
end
