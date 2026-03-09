import numpy as np
from Tools.Param import Parameters
from Dynamic.sysODE import cip, ssODEwrap, pend_cart_nl_wrap, wrapODE, pend_cart_nl, ssODE
from Tools.ctrl import lqrd, continuous_to_discrete
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from Nav.filter import ekf_update_d
from Guidance.trj import moveSidecip
from cip_acados import create_acados_ocp_cartpole
from acados_template import AcadosOcpSolver
from time import time

import os

print(os.getcwd())  # e.g., "C:/Users/you/project"
os.chdir("../cip_prj")
print(os.getcwd())  # Now prints: "C:/Users/you/project/cip_build"



ocp_solver = None


def u_lqr(x_cur, x_ref, Gain):
    """
    x_cur: current state vector (4,)
    x_ref:  guided state vector (4,)
    Gain:   error gain (1,4), in this case Klqr
    return: u = - Gain @ (x_curr - x_ref), (1,)
    """
    return -Gain @ (x_cur - x_ref)


def mpc_ctrl(x, ref, p):
    global ocp_solver

    dt = 0.01
    N = 20
    if ocp_solver is None:
        build_dir = "cip_build"
        dll_path = os.path.join(build_dir, "acados_ocp_solver_cartpole_tube_mpc.dll")
        json_path = 'acados_ocp.json'

        ocp = create_acados_ocp_cartpole(dt, N, p)
        ocp.code_export_directory = build_dir

        # === First time: generate code only ===
        if not os.path.exists(dll_path) or not os.path.exists(json_path):
            print("[INFO] Generating code in:", build_dir)
            AcadosOcpSolver(ocp, generate=True, build=False)
            print("[INFO] Code generated. Now run:")
            print(f"       mingw32-make -C {build_dir}")
            exit(1)

        # === After build: just load it ===
        ocp_solver = AcadosOcpSolver(ocp, json_file=json_path, generate=False, build=False)

    ocp_solver.set(0, "lbx", x)
    ocp_solver.set(0, "ubx", x)

    for i in range(ocp_solver.acados_ocp.dims.N):
        ocp_solver.set(i, "yref", np.concatenate([ref, [0.0]]))
    ocp_solver.set(N, "yref", ref)

    status = ocp_solver.solve()
    if status != 0:
        print(f"[WARN] Solver failed at step with status {status}")

    u0 = ocp_solver.get(0, "u")
    return u0


def show_results(tspan, x, u, perf):
    N = u.shape[0]
    nrows = x.shape[0] + 2
    fig, axes = plt.subplots(nrows=nrows, ncols=1, figsize=(8,8))
    plt.tight_layout(pad=3.)
    axes[0].plot(tspan, x[0], 'r+-', label='Position')
    axes[1].plot(tspan, x[1], 'b+-', label='Velocity')
    axes[2].plot(tspan, np.rad2deg(x[2]), 'r^-', label='Angle')
    axes[3].plot(tspan, np.rad2deg(x[3]), 'r^-', label='Rotation')
    axes[4].plot(tspan[:N], u, 'r^-', label='Control')
    axes[5].plot(tspan[:N], perf*1e3, 'b^-', label='Performance MPC')
    for ax in axes:
        ax.grid()
        ax.legend()
    plt.show(block=True)


def main():
    print("Checking Continious and Discrete system control")
    mu = 0.1962
    m = 0.248       # point mass
    M = 4.5236      # cart mass
    l = 0.765       # pendulum length
    b = 0.1         # linear damping
    d = 0.001       # rotational damping
    g = 9.80665     # g-Force
    J = mu * l ** 2 / 12
    mo = m + mu
    Mfull = M + mo
    Jc = J + mu * (l**2) / 4 + m * l**2
    Nt = 1000               # number of divisions
    dt = 0.01               # time step
    p = Parameters({        # make it possible to get fields via p.m, p.M, etc.
        "mu": mu,           # distributed mass
        "m": m,             # point mass
        "M": M,             # cart mass
        "Mfull": Mfull,     # cart mass
        "mo": mo,           # pendulum mass
        "l": l,             # pendulum length
        "b": b,             # linear damping
        "d": d,             # rotational damping
        "g": g,             # g-Force
        "J": J,
        "Jc": Jc,
        "MJ": Mfull * Jc - (mo * l) ** 2,
        "umin": -10,
        "umax": 10,
        "Nt": Nt,   # number of divisions
        "dt": dt,   # time step
        "w_max":    np.array([0.05, 0.1, np.deg2rad(2), np.deg2rad(3)]),
        "noise": 1,
    })
    pos0 = 0.; vel0 = 0.; eul0 = np.deg2rad(2.); rot0 = np.rad2deg(0.)
    x0 = np.array([pos0, vel0, eul0, rot0])
    u0 = np.array([0])
    nx = x0.shape[0]
    nu = u0.shape[0]
    ## LQR
    Ac, Bc = cip(x0, u0, p)                         # Continious system
    Ad, Bd = continuous_to_discrete(Ac, Bc, dt)     # Discrete system with dt time step
    Qc = np.eye(Ac.shape[1]) * 10
    Rc = np.array([[0.1]])
    K_lqrd, P = lqrd(Ac, Bc, Qc, Rc, dt)               # lqr control for discrete system
    ## Extended Kalman filter
    Ck = np.eye(nx)                                 # (4,4)
    Qk = np.diag([1e-4, 1e-3, 1e-4, 1e-3])          # Process noise covariance
    Rk = np.diag([1e-2, 1e-2, 1e-2, 1e-2])          # Measurement noise covariance
    P = np.eye(4) * 0.1

    # Initial condition
    tf = 10.
    tspan = np.linspace(0, tf, Nt + 1)
    tau = np.array([0, dt])
    x_real_hist = np.zeros((x0.shape[0], Nt+1))
    x_nom_hist = np.zeros_like(x_real_hist)
    x_est_hist = np.zeros_like(x_real_hist)
    x_ref_hist = np.zeros_like(x_real_hist)
    # Noise parameters
    w_max = np.array([0.05, 0.1, np.deg2rad(2), np.deg2rad(3)]) # noise amplitude for each channel (4, )
    w_a = 0.5
    # Storage
    x_est = x0 + w_a * w_max * (2 * np.random.randn(nx) - 1)
    x_real_hist[:, 0] = x0
    x_nom_hist[:, 0] = x0
    x_est_hist[:, 0] = x_est

    xr = np.zeros_like(x0)
    u_hist = np.zeros(Nt)
    u_prev = [0]
    tstop = 10.
    t_pause = 5.
    # Control constraints
    umin = -10
    umax = 10
    # Control rate constraints
    dumax = 100
    dumin = -100
    performance_mpc = np.zeros((Nt,))

    for i in range(Nt):
        # Get previous values
        xk = x_real_hist[:, i]
        xn = x_nom_hist[:, i]
        x_est = x_est_hist[:, i]

        # Guide
        # ref = moveSidecip(tspan[i], xr)         # Reference
        ref = xr
        start = time()
        # Define control
        u_nom = mpc_ctrl(xn, ref, p)       # MPC
        performance_mpc[i] = time() - start

        ui = u_lqr(x_est, xn, K_lqrd)  # LQR: -K_lqrd @ (x[:, i] - xr)

        ui += u_nom
        # Simulate plant behaviour
        res = odeint(lambda x, t: pend_cart_nl(x, ui[0], p), xk, tau)
        x_real_new = res[-1]

        # Update nominal system
        x_nom_new = Ad @ xn + Bd @ u_nom

        # Measurement with some noise
        y_meas = Ck @ x_real_new + w_max * np.random.multivariate_normal(np.zeros(nx), Rk)

        # Extended Kalman Filter update
        # xhat, P = ekf_update(xhat, P, ui, y_meas, dt, Qk, Rk, p)
        x_est, P = ekf_update_d(x_est, P, ui, y_meas, Ad, Bd, Ck, Qk, Rk)

        # Store results
        x_real_hist[:, i+1] = x_real_new
        x_nom_hist[:, i+1] = x_nom_new
        x_est_hist[:, i+1] = x_est
        x_ref_hist[:, i] = ref
        u_hist[i] = ui[0]
        u_prev = ui

    print("Solved")
    show_results(tspan, x_real_hist, u_hist, performance_mpc)


if __name__ == "__main__":
    main()