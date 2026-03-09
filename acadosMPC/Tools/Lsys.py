import numpy as np
from scipy.linalg import block_diag
from Tools.mathOp import Nw, Tw, grav_lin, skew
from Tools.tool import tensor_interpolation, mass_center_interpolate, ComponentForForces


def build_lsys(x, u, p):
    """
    Builds linearized system matrices A and B for a dynamical system.

    Parameters:
        x : ndarray - state vector for linearization
        u : ndarray - control vector
        p : object  - system parameters (must contain attributes:
                      m, mass, Inertia, geom, Cv, Cw, grav_on, g)

    Returns:
        A, B : tuple of ndarrays - system matrices
    """
    # Initialize matrices
    n_states = len(x)
    n_controls = len(u)
    A = np.zeros((n_states, n_states))
    B = np.zeros((n_states, n_controls))

    # Extract states
    r = x[0:3]  # coordinates
    v = x[3:6]  # velocity
    eul0 = x[6:9]  # Euler angles
    w = x[9:12]  # angular velocity

    # Unpack parameters
    m = p.m
    Mass = np.diag([m, m, m])
    Js = tensor_interpolation(m, p.mass, p.Inertia)
    rcg = mass_center_interpolate(m, p.mass, p.rcgs)
    Hf, Hm = ComponentForForces(p.geom, rcg)
    Cv = p.Cv
    Cw = p.Cw

    # Build A matrix blocks
    A11 = np.zeros((3, 3))
    A12 = np.eye(3)
    A21 = np.zeros((3, 3))
    A22 = -skew(w) + np.linalg.solve(Mass, -Cv)

    At = np.block([[A11, A12],
                   [A21, A22]])

    A33 = np.zeros((3, 3))
    A43 = np.zeros((3, 3))
    A34 = Tw(eul0)  # Assuming tw() is equivalent to MATLAB's Tw()
    A44 = np.linalg.solve(Js, -skew(w) @ Js - Cw)

    Ar = np.block([[A33, A34],
                   [A43, A44]])

    # Combine A matrix
    A = block_diag(At, Ar)

    # Add gravity terms if enabled
    if p.grav_on > 0:
        A[3:6, 6:9] = grav_lin(eul0, p.g)

    # Build B matrix
    B11 = np.zeros((3, 9))
    B21 = np.linalg.solve(Mass, Hf)
    Bt = np.vstack([B11, B21])

    B31 = np.zeros((3, 9))
    B41 = np.linalg.solve(Js, Hm)
    Br = np.vstack([B31, B41])

    B = np.vstack([Bt, Br])

    return A, B

