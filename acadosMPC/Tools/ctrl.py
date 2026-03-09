import numpy as np
from scipy.optimize import minimize
import scipy.linalg as la
from scipy.linalg import expm


# from scipy.signal import cont2discrete
# A_d, B_d, _, _ = cont2discrete((A, B, np.eye(2), np.zeros((4,1))), dt, method='zoh')


def lqr(a, b, q, r):
    """
    xdot = Ax + Bu
    x - state vector (N x 1)
    a - system matrix (N x 1)
    b - control matrix (N x M)
    q - weight errors coeffs matrix (N x N)
    r - weight control coeff matrix (M x 1)
    """
    p = la.solve_continuous_are(a, b, q, r)
    # k =  np.linalg.inv(r) * b @ p - в общем виде
    k = la.inv(r) * b.T @ p             # в случае если r - число (или матрица из одного элемента)
    return k, p


def dlqr(a, b, q, r):
    """
    Discrete-time Linear Quadratic Regulator

    x[k+1] = A x[k] + B u[k]
    x - state vector (n x 1)
    a - system matrix (n x n)
    b - control matrix (n x m)
    q - state cost matrix (n x n)
    r - control cost matrix (m x m)
    """
    # Solve the discrete algebraic Riccati equation
    p = la.solve_discrete_are(a, b, q, r)

    # Compute the LQR gain
    k = la.inv(b.T @ p @ b + r) @ (b.T @ p @ a)

    return k, p


def lqrd(A, B, Q, R, dt=None):
    """
    Discrete-time LQR controller.

    Args:
        A: State matrix (n x n) — discrete or continuous.
        B: Control matrix (n x m) — discrete or continuous.
        Q: State cost matrix (n x n).
        R: Control cost matrix (m x m).
        dt: Time step (if A, B are continuous-time).
    Returns:
        K: Optimal gain matrix (m x n).
    """
    if dt is not None:
        A, B = continuous_to_discrete(A, B, dt)
        Q *= dt
        R *= dt

    # Solve Discrete Algebraic Riccati Equation (DARE)
    P = la.solve_discrete_are(A, B, Q, R)

    # Compute optimal gain K
    K = la.inv(B.T @ P @ B + R) @ (B.T @ P @ A)
    return K, P


def continuous_to_discrete(A_c, B_c, dt):
    n = A_c.shape[0]
    M = np.block([[A_c, B_c], [np.zeros((B_c.shape[1], n + B_c.shape[1]))]])
    expM = expm(M * dt)
    A_d = expM[:n, :n]
    B_d = expM[:n, n:]
    return A_d, B_d


def approximate_rpi_box(A_K, nx, w_bound, max_iter=100, tol=1e-6):
    """
    Approximates the Robust Positively Invariant (RPI) set as a box.

    Args:
        A_K: Closed-loop system matrix (A - BK) (nx x nx).
        nx: Number of states.
        w_bound: Disturbance bound (scalar or vector).
        max_iter: Maximum iterations (default: 100).
        tol: Convergence tolerance (default: 1e-6).

    Returns:
        Z_box: RPI box bounds (nx x 1 vector).
    """
    Z = np.eye(nx)
    for k in range(max_iter):
        Z_new = A_K @ Z + np.eye(nx)  # Unit box approximation for disturbance
        if np.max(np.abs(Z_new - Z)) < tol:
            break
        Z = Z_new
    Z_box = w_bound * np.max(np.abs(Z), axis=1)  # Conservative bounding box
    return Z_box


def pid_controller(xh, xr, p, integral_sum=None, prev_error=None, mode=0):
    """
    PID Controller implementation for 3D vector control

    Parameters:
    - xh: current value (3,) array
    - xr: reference vector (3,) array
    - p: parameter object with p.ctrl.p, p.ctrl.i, p.ctrl.d as (3,1) arrays
    - dt: time step for derivative and integral calculations
    - integral_sum: accumulated integral error (3,1), None for first call
    - prev_error: previous error for derivative calculation (3,1), None for first call

    Returns:
    - control_vector: PID control output (3,)
    - current_error: current error (3,)
    - integral_sum: updated integral sum for next iteration
    """

    # Ensure proper shapes
    # xh = np.array(xh).reshape(3, )  # Convert (3,) to (3,1)
    # xr = np.array(xr).reshape(3, 1)  # Ensure (3,1) shape
    # Calculate current error
    current_error = xr - xh

    dt = p.ctrl.dt
    # Proportional term
    proportional = p.ctrl.p * current_error

    # Integral term
    if integral_sum is None:
        integral_sum = np.zeros(3)
    integral_sum += current_error * dt
    integral = p.ctrl.i * integral_sum

    # Derivative term
    if prev_error is None:
        derivative = np.zeros(3)
    else:
        derivative = p.ctrl.d * (current_error - prev_error) / dt

    # Calculate control vector
    control_vector = proportional + integral + derivative
    if mode == 0:
        u = np.hstack((control_vector, np.zeros(3)))
    elif mode == 1:
        u = np.hstack((np.zeros(3), control_vector))
    else:
        u = control_vector
    return u, current_error, integral_sum



def RPI_Check():
    # Define system matrices (example: double integrator)
    Ad = np.array([[1, 0.1], [0, 1]])
    Bd = np.array([[0], [1]])
    K_lqr = np.array([[0.5, 1.0]])  # LQR gain (1 x nx)

    # Closed-loop matrix
    A_K = Ad - Bd @ K_lqr

    # Parameters
    nx = 2
    w_max = 0.01  # Disturbance bound

    # Compute RPI box
    Z_rpi_bounds = approximate_rpi_box(A_K, nx, w_max, max_iter=100, tol=1e-6)
    print("RPI box bounds:", Z_rpi_bounds)
    delta_u = abs(K_lqr @ Z_rpi_bounds)
    print("delta_u:", delta_u)



if __name__ == "__main__":
    A_c = np.array([[0, 1], [0, 0]])  # Example: Double integrator
    B_c = np.array([[0], [1]])
    dt = 0.2
    # A_d, B_d, _, _, dt1 = cont2discrete((A_c, B_c, np.eye(2), np.zeros((2,1))), dt, method='zoh')
    A2, B2 = continuous_to_discrete(A_c, B_c, dt)
    print("A")
    print(A2)
    print("B")
    print(B2)
    RPI_Check()