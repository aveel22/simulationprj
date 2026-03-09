import numpy as np


def cascade_ODE(t, xs, u, T, n):
    dx = np.zeros_like(xs)
    if xs[-1] > u:
        T /= 4.
    dx[0] = (u - xs[0]) / T
    for i in range(1, n):
        dx[i] = (xs[i - 1] - xs[i]) / T
    return dx


def transition_ODE(t, xs, u, L, n):
    """
    L: Rising Time
    n: ODE's Order
    """
    return cascade_ODE(t, xs, u, L / n, n)
