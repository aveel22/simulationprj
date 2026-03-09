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
import subprocess
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
    N = 30
    if ocp_solver is None:
        build_dir = "cip_build"
        dll_path = os.path.join(build_dir, "acados_ocp_solver_cartpole_mpc.dll")
        json_path = 'acados_ocp.json'

        ocp = create_acados_ocp_cartpole(dt, N, p)
        ocp.code_export_directory = build_dir

        # === First time: generate code only ===
        if not os.path.exists(dll_path) or not os.path.exists(json_path):
            print("[INFO] Generating code in:", build_dir)
            AcadosOcpSolver(ocp, json_file=json_path, generate=True, build=False)
            print("[INFO] Code generated. Now run:")
            print(f"       mingw32-make -C {build_dir}")
            exit(1)

        # === After build: just load it ===
        ocp_solver = AcadosOcpSolver(ocp, json_file=json_path, generate=False, build=False)
        # subprocess.check_call(["mingw32-make", "-C", build_dir])

    ocp_solver.set(0, "lbx", x.copy())
    ocp_solver.set(0, "ubx", x.copy())

    # Optional warm start
    if ocp_solver.acados_ocp.dims.N > 1:
        ocp_solver.set(1, "x", x.copy())  # ← warm start guess for the second shooting node

    for i in range(ocp_solver.acados_ocp.dims.N):
        ocp_solver.set(i, "yref", np.concatenate([ref, [0.0]]))

    ocp_solver.set(N, "yref", ref.copy())  # ✅ terminal reference

    status = ocp_solver.solve()
    if status != 0:
        print(f"[WARN] Solver failed at step with status {status}")

    u0 = ocp_solver.get(0, "u")
    return u0


def show_results(tspan, x, options=None, title=None):

    nrows = x.shape[0]
    fig, axes = plt.subplots(nrows=nrows, ncols=1, figsize=(8,8))
    plt.tight_layout(pad=3.)
    if title:
        plt.title(title)
    if nrows == 1:
        if options:
            axes.plot(tspan, x[0], 'r+-', label=options[0]["label"])
            axes.set_xlabel(options[0]["xlabel"])
            axes.set_ylabel(options[0]["ylabel"])
        else:
            axes.plot(tspan, x[0], 'r+-', label=f'arg #1')
        axes.grid()
        axes.legend()
    else:
        for i, ax in enumerate(axes):
            if options:
                ax.plot(tspan, x[i], 'r+-', label=options[i]["label"])
                ax.set_xlabel(options[i]["xlabel"])
                ax.set_ylabel(options[i]["ylabel"])
            else:
                ax.plot(tspan, x[i], 'r+-', label=f"# arg{i+1}")
            ax.legend()
            ax.grid()



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
        "umin": -20,
        "umax": 20,
        "dumin": -500,
        "dumax": 500,
        "Nt": Nt,   # number of divisions
        "dt": dt,   # time step
        "w_max":    np.array([0.05, 0.1, np.deg2rad(0.5), np.deg2rad(0.5)]),
        "noise": 0,
    })
    pos0 = 0.; vel0 = 0.; eul0 = np.deg2rad(2.); rot0 = np.rad2deg(0.)
    x0 = np.array([pos0, vel0, eul0, rot0])
    u0 = np.array([0])
    nx = x0.shape[0]
    nu = u0.shape[0]
    ## LQR
    Ac, Bc = cip(x0, u0, p)                         # Continious system
    Ad, Bd = continuous_to_discrete(Ac, Bc, dt)     # Discrete system with dt time step
    Qc = np.diag([100., 10., 100., 10.])
    Rc = np.array([[0.01]])
    K_lqrd, Plqr = lqrd(Ac, Bc, Qc, Rc, dt)               # lqr control for discrete system
    ## Extended Kalman filter
    Ck = np.eye(nx)                                 # (4,4)
    Qk = np.diag([1e-4, 1e-3, 1e-4, 1e-3])          # Process noise covariance
    Rk = np.diag([1e-2, 1e-2, 1e-2, 1e-2])          # Measurement noise covariance
    P = np.eye(4) * 0.1
    p['Ad'] = Ad
    p['Bd'] = Bd
    p['P'] = Plqr
    # Initial condition
    tf = Nt * dt
    tspan = np.linspace(0., tf, Nt+1)
    tau = np.array([0, dt])
    xk_hist = np.zeros((x0.shape[0], Nt+1))
    xn_hist = np.zeros_like(xk_hist)
    x_est_hist = np.zeros_like(xk_hist)
    x_ref_hist = np.zeros_like(xk_hist)
    # Noise parameters
    w_max = np.array([0.05, 0.1, np.deg2rad(0.5), np.deg2rad(0.5)]) # noise amplitude for each channel (4, )
    w_a = 0.5
    # Storage
    x_est = x0 + w_a * w_max * (2 * np.random.randn(nx) - 1)
    xk_hist[:, 0] = x0
    xn_hist[:, 0] = x0
    x_est_hist[:, 0] = x_est

    xr = np.zeros_like(x0)
    u_hist = np.zeros(Nt+1)
    performance_mpc = np.zeros((Nt+1,))

    for i in range(Nt):
        # Get previous values
        xk = xk_hist[:, i]
        # xn = x_nom_hist[:, i]

        # Guide
        ref = xr
        # Define control
        start = time()
        u_ = mpc_ctrl(xk, ref, p)       # MPC
        performance = time() - start
        ui = u_

        # Simulate plant behaviour
        res = odeint(lambda x, t: pend_cart_nl(x, ui[0], p), xk, tau)
        xk_new = res[-1]

        # Store results
        xk_hist[:, i+1] = xk_new
        x_ref_hist[:, i] = ref
        performance_mpc[i] = performance
        u_hist[i] = ui[0]


    print("Solved")
    xk_hist[2:] = np.rad2deg(xk_hist[2:])
    opt_1 = {0: {"xlabel": "Time, s", "ylabel": "X [m]", "label": "Position"},
             1: {"xlabel": "Time, s", "ylabel": "V [m/s]", "label": "Velocity"},
             2: {"xlabel": "Time, s", "ylabel": "$\Theta$"+" [deg]", "label": "Angle"},
             3: {"xlabel": "Time, s", "ylabel": "$\omega$"+" [deg/s]", "label": "Rotation"},
             }
    opt_2 = {0: {"xlabel": "Time, s", "ylabel": "Perf. [ms]", "label": "Performance"}}
    opt_3 = {0: {"xlabel": "Time, s", "ylabel": "F [N]", "label": "Force"}}
    show_results(tspan, np.vstack((xk_hist, )), options=opt_1, title="State")
    show_results(tspan, np.vstack((1e3*performance_mpc, )), options=opt_2, title="Calculation Performance")
    show_results(tspan, np.vstack((u_hist, )), options=opt_3, title="Control")
    plt.show(block=True)


if __name__ == "__main__":
    main()
