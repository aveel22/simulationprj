import numpy as np
from Dynamic.sysODE import cip


def ekf_update(xhat, P, u, y, dt, Q, R, p):
    """
    Extended Kalman Filter (EKF) update step for a nonlinear system.

    Args:
        xhat: Current state estimate (4x1 vector).
        P: Current error covariance matrix (4x4).
        u: Control input (scalar or 1x1).
        y: Measurement vector (4x1, assuming full-state feedback).
        dt: Time step.
        Q: Process noise covariance (4x4).
        R: Measurement noise covariance (4x4).
        p: System parameters (passed to pendcart_nl and estimate_A).

    Returns:
        xhat_new: Updated state estimate.
        P_new: Updated error covariance.
    """
    # --- Prediction Step ---
    # Compute nonlinear dynamics: xdot = f(x, u)

    # Linearize A matrix (Jacobian of f(x, u) w.r.t. x)
    A, B = cip(xhat,u,p)

    # Predict state (Euler integration)
    xhat_pred = xhat + A @ xhat * dt

    # Predict covariance
    P_pred = A @ P @ A.T + Q

    # --- Update Step ---
    # Measurement matrix (identity since y = x + noise)
    C = np.eye(4)

    # Kalman gain
    S = C @ P_pred @ C.T + R
    Kf = P_pred @ C.T @ np.linalg.inv(S)

    # Update state estimate
    y_residual = y - C @ xhat_pred
    xhat_new = xhat_pred + Kf @ y_residual

    # Update covariance (Joseph form for stability)
    I = np.eye(4)
    P_new = (I - Kf @ C) @ P_pred @ (I - Kf @ C).T + Kf @ R @ Kf.T

    return xhat_new, P_new


def ekf_update_d(xhat, P, u, y, Ad, Bd, C_meas, Q_proc, R_meas):
    """
    Discrete Extended Kalman Filter (EKF) update step.

    Args:
        xhat: Current state estimate (nx,).
        P: Current error covariance (nx, nx).
        u: Control input (nu,).
        y: Measurement vector (ny,).
        Ad: Discrete state transition matrix (nx, nx).
        Bd: Discrete control input matrix (nx, nu).
        C_meas: Measurement matrix (ny, nx).
        Q_proc: Process noise covariance (nx, nx).
        R_meas: Measurement noise covariance (ny, ny).

    Returns:
        xhat_new: Updated state estimate.
        P_new: Updated error covariance.
    """
    # --- Prediction Step ---
    x_pred = Ad @ xhat + Bd @ u
    P_pred = Ad @ P @ Ad.T + Q_proc

    # --- Update Step ---
    S = C_meas @ P_pred @ C_meas.T + R_meas
    K_kalman = P_pred @ C_meas.T @ np.linalg.inv(S)
    y_residual = y - C_meas @ x_pred

    xhat_new = x_pred + K_kalman @ y_residual
    P_new = (np.eye(len(xhat)) - K_kalman @ C_meas) @ P_pred

    return xhat_new, P_new