import numpy as np



def rotx(a):
    # % ROTX Build rotation matrix around X axis
    # % x - vector
    return np.array([
        [1,     0,      0],
        [0, np.cos(a), np.sin(a)],
        [0, -np.sin(a), np.cos(a)],
    ])


def roty(a):
    # % ROTX Build rotation matrix around X axis
    # % x - vector
    return np.array([
        [np.cos(a), 0, -np.sin(a)],
        [0, 1, 0],
        [np.sin(a), 0, np.cos(a)],
    ])


def rotz(a):
    # % ROTX Build rotation matrix around X axis
    # % x - vector
    return np.array([
        [np.cos(a), np.sin(a), 0],
        [-np.sin(a), np.cos(a), 0],
        [0, 0, 1],
    ])

# --- derivatives consistent with the above definitions ---
def drotx_da(a):
    c = np.cos(a); s = np.sin(a)
    # d/da of rotx:
    # [[0,0,0],
    #  [0,-sin, cos],
    #  [0,-cos,-sin]]
    return np.array([
        [0, 0, 0],
        [0, -s,  c],
        [0, -c, -s],
    ])

def droty_da(a):
    c = np.cos(a); s = np.sin(a)
    # d/da of roty:
    # [[-sin, 0, -cos],
    #  [    0, 0,     0],
    #  [ cos, 0, -sin]]
    return np.array([
        [-s, 0, -c],
        [ 0, 0,  0],
        [ c, 0, -s],
    ])

def drotz_da(a):
    c = np.cos(a); s = np.sin(a)
    # d/da of rotz:
    # [[-sin, cos, 0],
    #  [-cos,-sin, 0],
    #  [   0,   0, 0]]
    return np.array([
        [-s,  c, 0],
        [-c, -s, 0],
        [ 0,  0, 0],
    ])



def Nw(eul):
    # This is transformation matrix from Launch to Body Frame
    # %Pitch(Z)->Yaw(Y)->Roll(X)
    # Unpack the input vector
    roll = eul[0]
    pitch = eul[1]
    yaw = eul[2]
    return rotx(roll) @ roty(yaw) @ rotz(pitch)


def skew(x):
    return np.array([
        [0, -x[2], x[1]],
        [x[2], 0, -x[0]],
        [-x[1], x[0], 0]
    ])


def Tw(eul):
    # Unpack the input vector
    # Rocket transformation
    # Pitch(Z) -> Yaw(Y) -> Roll(Z)
    # Use of matrix
    # [phi; theta; psi] = Tw(eul) * [wx; wy; wz]

    roll = eul[0]        # roll
    pitch = eul[1]      # pitch
    yaw = eul[2]        # yaw

    # Check for singularity(pitch=pi / 2 or -pi / 2 leads to division by zero)
    if np.cos(yaw) == 0:
        raise Exception('pitch must not be +/- pi/2 to avoid division by zero.')

    # Compute the transformation matrix
    return np.array([
        [np.cos(yaw), np.sin(yaw) * np.sin(roll), np.sin(yaw) * np.cos(roll)],
        [0, np.sin(roll), np.cos(roll)],
        [0, np.cos(yaw) * np.cos(roll), -np.cos(yaw) * np.sin(roll)]
    ]) / np.cos(yaw)


def E2W(eul):
    # E2W Summary of this function goes here
    # Euler's kinematic for NED
    # Yaw(Z) -> Pitch(Y) -> Roll(Z)
    # eul:      roll, pitch, yaw [radians]
    # returns:  K -> eul_dot = K @ [wx, wy, wz]
    # roll = eul[0]
    # pitch = eul[1]
    # yaw = eul[2]
    # return np.array([
    #         [np.cos(pitch), np.sin(roll) * np.sin(pitch), np.cos(roll) * np.sin(pitch)],
    #         [0, np.cos(roll) * np.cos(pitch), -np.sin(roll) * np.cos(pitch)],
    #         [0, np.sin(roll), np.cos(roll)]
    #     ] ) / np.cos(pitch)

    s = np.sin(eul)
    c = np.cos(eul)
    return np.array([
            [c[1], s[0] * s[1],  c[0] * s[1]],
            [0,    c[0] * c[1], -s[0] * c[1]],
            [0,           s[0],         c[0]]
        ] ) / c[1]


def K2W(eul):
    # E2W Summary of this function goes here
    # Euler's kinematic for ENU
    # Yaw(Z) -> Pitch(Y) -> Roll(Z)
    # Unlike NED it has opposite pitch rotation
    return E2W([eul[0], -eul[1], eul[2]])


def W2E(eul):
    """
    Euler's kinematic for NED CS
    eul:    roll, pitch, yaw [radians]
    return: kinematic matrix K, e.g. w = K @ eul_dot
    """
    s = np.sin(eul)
    c = np.cos(eul)
    return np.array([
        [1,     0,        -s[1]],
        [0,  c[0],  s[0] * c[1]],
        [0, -s[0],  c[0] * c[1]],
    ])

def W2K(eul):
    # E2W Summary of this function goes here
    # Euler's kinematic for ENU
    # Yaw(Z) -> Pitch(Y) -> Roll(Z)
    # Unlike NED it has opposite pitch rotation
    return W2E(np.array([eul[0], -eul[1], eul[2]]))

def grav_lin(eul, g_c):
    # Возвращает 3x3 матрицу производных ускорения от гравитации
    # по углам arg = [phi, theta, psi] в body frame
    # Параметры:
    #   phi   — угол крена (в радианах)
    #   theta — угол тангажа (в радианах)
    #   psi   — угол рыскания (в радианах)
    #   g_c   — ускорение свободного падения (>0) in body frame
    #
    # Выход:
    #   A     — 3x3 матрица: d(ag_body)/d(phi, theta, psi)
    #

    phi = eul[0]        #  roll
    theta = eul[1]      # pitch
    psi = eul[2]        # yaw
    g = g_c


    # % Предварительные значения
    sphi = np.sin(phi);   cphi = np.cos(phi)
    stheta = np.sin(theta); ctheta = np.cos(theta)
    spsi = np.sin(psi);   cpsi = np.cos(psi)

    # Строки по производным
    A11 = 0
    A12 = -g * cpsi * ctheta
    A13 =  g * spsi * stheta

    A21 =  g * (sphi * ctheta - spsi * stheta * cphi)
    A22 =  g * (-sphi * spsi * ctheta + stheta * cphi)
    A23 = -g * sphi * stheta * cpsi

    A31 =  g * (sphi * spsi * stheta + cphi * ctheta)
    A32 = -g * (sphi * stheta + spsi * cphi * ctheta)
    A33 = -g * stheta * cphi * cpsi

    # % Сборка в матрицу
    return np.array([
        [A11, A12, A13],
        [A21, A22, A23],
        [ A31, A32, A33]
    ])

def N2B(eul):
    #  This is transformation matrix from Launch to Body Frame
    #  Yaw(Z)->Pitch(Y)->Roll(X)
    # Unpack the input vector
    roll = eul[0]       # roll
    pitch = -eul[1]     # pitch
    yaw = eul[2]        # yaw
    Nrot = rotx(roll) @ roty(pitch) @ rotz(yaw)
    return Nrot


def D2B_from_euler(eul, g):
    """
    Returns the derivative of N2B(eul) w.r.t. Euler angles,
    so that (D2B_from_euler(eul, g) = linearized gravity acceleration.
    """
    roll, pitch, yaw = eul

    Rx = rotx(roll); Ry = roty(pitch); Rz = rotz(yaw)
    dRx = drotx_da(roll); dRy = droty_da(pitch); dRz = drotz_da(yaw)

    # derivative contributions
    dR_droll = dRx @ Ry @ Rz
    # ∂R/∂pitch = Rx * dRy * Rz
    dR_dpitch = Rx @ dRy @ Rz
    # ∂R/∂yaw   = Rx * Ry * dRz
    dR_dyaw = Rx @ Ry @ dRz

    c_roll = dR_droll @ g
    c_pitch = dR_dpitch @ g
    c_yaw = dR_dyaw @ g

    D2B = np.column_stack([c_roll, c_pitch, c_yaw])
    # Pack as 3x3x3 tensor: each slice is ∂N/∂eul_i
    return D2B


def D2B_numerical(eul, g, eps=1e-6):
    # finite-difference check
    D_num = np.zeros((3, 3))
    g0 = N2B(eul) @ g

    for i in range(3):
        de = np.zeros(3)
        de[i] = eps
        g1 = N2B(eul + de) @ g
        D_num[:, i] = (g1 - g0) / eps

    return D_num

def check_gravity_linearization(run=False):
    if not run:
        return
    eul1 = np.deg2rad([5.0, 2.0, 10.0])         # roll, pitch, yaw
    eul = np.deg2rad([0.0, 0.0, 0.0])           # roll, pitch, yaw
    g_nav = np.array([-9.80665, 0.0, 0])      # gravity in nav frame

    # analytic D2B (your convention)
    D_analytic = D2B_from_euler(eul, g_nav)
    print("D_analytic:\n", D_analytic)

    # finite-difference check
    eps = 1e-6
    D_num = np.zeros((3, 3))
    R0 = N2B(eul)
    g0 = R0 @ g_nav

    for i in range(3):
        de = np.zeros(3)
        de[i] = eps
        R1 = N2B(eul + de)
        g1 = R1 @ g_nav
        D_num[:, i] = (g1 - g0) / eps

    print("\nD_numeric (finite diff):\n", D_num)
    print("\nDifference (analytic - numeric):\n", D_analytic - D_num)


def check_N2B(run=False):
    if not run:
        return
    eul = np.array([0, np.pi/2, np.pi/2])
    eul1 = np.array([0, 0, np.pi/2])
    eul2 = np.array([0, np.pi/2, np.pi/2])
    print(f"Result 0: {eul}")
    print(Nw(eul), end="\n\n")
    print(N2B(eul))
    print(f"Result 1: {eul1}")
    print(Nw(eul1), end="\n\n")
    print(N2B(eul1))
    print(f"Result 2: {eul2}")
    print(Nw(eul2))
    print(N2B(eul2))


def check_K2W(run=False):
    w = np.array([0, 10, 0])
    eul = np.deg2rad([0, 90, 0])
    eul_dot = np.array([10, 10, 10])
    eul_dot2 = np.array([10, 10, -10])
    # print(K2W(eul) @ w)
    # print(E2W(eul) @ w)
    print("W2E = \n", W2E(eul))
    print("W2K = \n", W2K(eul))
    print("*"*20)
    print(W2E(eul) @ eul_dot)
    print(W2K(eul) @ eul_dot2)

    print("*=*" * 20)
    from scipy.integrate import solve_ivp
    import scipy.linalg as la
    import matplotlib.pyplot as plt

    def ode(x, u, p):
        J = p["J"]
        Cw = p["Cw"]
        return np.hstack([
            K2W(x[:3]) @ x[3:],
            la.inv(J) @ (-Cw @ x[3:] + u)
        ])

    def ode2(x, u, p):
        J = p["J"]
        Cw = p["Cw"]
        return np.hstack([
            E2W(x[:3]) @ x[3:],
            la.inv(J) @ (-Cw @ x[3:] + u)
        ])


    ndim = eul.shape[0]
    Cw = np.eye(ndim) * 0.1
    T0, Tf, dt = 0., 1, 0.1
    tspan = np.arange(T0, Tf+dt, dt)
    Nt = tspan.shape[0]

    eul = np.zeros((int(2 * ndim), Nt))
    eul2 = np.zeros_like(eul)
    J = np.eye(ndim) * 10.
    M = np.array([2, 1, 10])
    p = {"J": J, "Cw": Cw}
    for i in range(Nt-1):
        tau = tspan[i:i+2]
        u = M
        sol = solve_ivp(lambda t, x : ode(x,u,p), tau, eul[:, i], t_eval=[tau[1]])
        sol2 = solve_ivp(lambda t, x : ode2(x,u,p), tau, eul2[:, i], t_eval=[tau[1]])
        eul[:, i+1] = sol.y[:, -1]
        eul2[:, i+1] = sol2.y[:, -1]




    labels = ["$\\varphi$", "$\\vartheta$", "$\psi$", "$\omega_X$", "$\omega_Y$", "$\omega_Z$"]
    yaxis = ["Angles [DEG]", "Rates [DEG/S]"] #, "$Moment [N \cdot m]$"]
    fig, axes = plt.subplots(nrows=2,ncols=2, figsize=(14, 8))
    eul = np.rad2deg(eul)
    eul2 = np.rad2deg(eul2)
    for i, row in enumerate(axes):
        for j in range(3):
            row[0].plot(tspan, eul[3*i+j], label=labels[3*i+j])
            row[1].plot(tspan, eul2[3*i+j], label=labels[3*i+j])
        for ax in row:
            ax.set_xlabel("Time [S]")
            ax.set_ylabel(yaxis[i])
            ax.legend()
            ax.grid()

    plt.show()


if __name__ == "__main__":
    check_N2B(run=False)
    check_gravity_linearization(run=True)
    check_K2W(run=False)