import numpy as np
from Tools.Param import Parameters
import casadi as ca


class ODEAdapter:
    def __init__(self, func):
        self.func = func
        self.call_count = 0  # Can track stats

    def __call__(self, x, t, *args):
        self.call_count += 1
        u, p = args
        return self.func(x, u, p)


# Apply as decorator
@ODEAdapter
def pend_cart_nl2(x, u, p):
    return pend_cart_nl(x, u, p)

"""
# Usage
ode_func = pend_cart_nl  # Now odeint-compatible
sol = odeint(ode_func, x0, t, args=(u, p))
print(f"ODE was called {ode_func.call_count} times")
"""

class wrapODE:
    """
    Too complicated
    """
    def __call__(self, func, q, t, *args):
        return func(q, t, *args)


def ssODEwrap(q, t, u, A, B):
    x = q
    return ssODE(x, u, A, B)


def ssODE(x, u, A, B):
    """
    ODE for Linear system
    """
    return A @ x + B @ u


def cip(x,u,p):
    """
    Cart-Pole Linear system. Return A, B matrices
    m = 0.248  # point mass
    M = 4.5236  # cart mass
    l = 0.765  # pendulum length
    b = 0.5  # linear damping
    d = 0.018  # rotational damping
    g = 9.80665  # g-Force
    Nt = 1000  # number of divisions
    dt = 0.01  # time step
    """
    b = p.b
    Jc = p.Jc
    MJ = p.MJ
    mo = p.mo
    L = p.l
    g = p.g
    Mfull = p.Mfull
    Nt = p.Nt           # number of divisions
    dt = p.dt           # time step

    A = np.array([
        [0, 1, 0, 0],
        [0, -b*Jc / MJ, (mo*L) ** 2 * g / MJ, 0],
        [0, 0, 0, 1],
        [0, -b*mo*L / MJ, Mfull*mo*L*g / MJ, 0]
    ])

    B = np.array([[0], [Jc / MJ], [0], [mo*L / MJ]])  # Control Matrix
    return A, B


def cip_variation(x, u, p, f):
    new_params = Parameters({
        "b": p.b * (1 + f * (2 * np.random.rand() - 1)),
        "Jc": p.Jc * (1 + f * (2 * np.random.rand() - 1)),
        "MJ": p.MJ * (1 + f * (2 * np.random.rand() - 1)),
        "mo": p.mo * (1 + f * (2 * np.random.rand() - 1)),
        "l": p.l * (1 + f * (2 * np.random.rand() - 1)),
        "g": p.g * (1 + f * (2 * np.random.rand() - 1)),
        "Mfull": p.Mfull * (1 + f * (2 * np.random.rand() - 1)),
    })
    return cip(x, u, new_params)


def pend_cart_nl_wrap(x, t, u, p):
    w = p.w_max * (2 * np.random.randn(4) - 1) * p.noise
    return pend_cart_nl(x, u[0], p) + w


def pend_cart_nl(x, u, p):
    """
    x[0]    - x
    x[1]    - v
    x[2]    - phi
    x[3]    - omega
    u       - F
    """
    # Nonlinear Dynamics function
    g = p.g
    b = p.b
    M = p.Mfull
    mo = p.mo
    Jc = p.Jc
    L = p.l

    Sx = np.sin(x[2])
    Cx = np.cos(x[2])
    D = M * Jc - (mo * L * Cx) ** 2

    return np.array([
        x[1],
        (1 / D) * (-(b * x[1] + mo * L * Sx * x[3] ** 2) * Jc + (mo * L) ** 2 * g * Sx * Cx + Jc * u),
        x[3],
        (1 / D) * (-mo * L * Cx * (b * x[1] + mo * L * Sx * x[3] ** 2) + M * mo * g * L * Sx + mo * L * Cx * u)
    ])


def pend_cart_nl_sym(x, u, p):
    """
    Symbolic dynamics for the cart-pole system, for use in acados.
    Arguments:
        x: casadi SX (4x1) - [x, v, phi, omega]
        u: casadi SX (1x1) - force input
        p: Parameters class with fields: g, b, Mfull, mo, Jc, l
    Returns:
        x_dot: casadi SX (4x1)
    """

    # Unpack symbolic states
    pos = x[0]
    vel = x[1]
    phi = x[2]
    omega = x[3]

    # Trig functions
    Sx = ca.sin(phi)
    Cx = ca.cos(phi)

    # Parameters
    g = p.g
    b = p.b
    M = p.Mfull
    mo = p.mo
    Jc = p.Jc
    L = p.l  # Make sure it's lowercase "l" if that's how you store it

    # Denominator
    D = M * Jc - (mo * L * Cx)**2

    # Dynamics
    x_dot = ca.vertcat(
        vel,
        (1 / D) * (-(b * vel + mo * L * Sx * omega**2) * Jc + (mo * L)**2 * g * Sx * Cx + Jc * u),
        omega,
        (1 / D) * (-mo * L * Cx * (b * vel + mo * L * Sx * omega**2) + M * mo * g * L * Sx + mo * L * Cx * u)
    )
    return x_dot


"""
x_sym = ca.SX.sym('x', 4)
u_sym = ca.SX.sym('u', 1)
x_dot_sym = pend_cart_nl_sym(x_sym, u_sym, p)
dyn_fun = ca.Function('pend_cart_nl_fun', [x_sym, u_sym], [x_dot_sym])
"""