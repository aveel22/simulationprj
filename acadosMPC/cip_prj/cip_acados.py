from acados_template import AcadosOcp, AcadosModel
from casadi import SX, vertcat
from Dynamic.sysODE import pend_cart_nl_sym
import numpy as np


def create_acados_ocp_cartpole(dt, N, p):
    model = AcadosModel()
    model.name = 'cartpole_mpc'

    # States
    x = SX.sym('x')
    v = SX.sym('v')
    th = SX.sym('th')
    w = SX.sym('w')
    x_ = vertcat(x, v, th, w)

    # Control
    u = SX.sym('u')

    # Cost reference dimensions
    nu = 1
    nx = 4
    ny = nx + nu  # 4 states + 1 input

    # Dynamics (you can use symbolic version or cip() wrapper)
    f = pend_cart_nl_sym(x_, u, p)  # You can use casadi model directly or rewrap

    model.x = x_
    model.u = u
    model.f_expl_expr = f
    model.f_impl_expr = x_ - f  # for implicit

    ocp = AcadosOcp()
    ocp.model = model

    # Horizon and time step
    # ocp.dims.N = N
    # ocp.solver_options.N_horizon = N
    ocp.dims.N = N  # ✅ CORRECTED — needed to define N here explicitly
    ocp.solver_options.tf = dt * N

    # Cost
    Q = np.diag([100, 10, 100, 10])
    R = np.array([[0.01]])
    W = np.block([
        [Q, np.zeros((4, 1))],
        [np.zeros((1, 4)), R]
    ])
    # === Cost setup ===
    ocp.cost.cost_type_0 = 'LINEAR_LS'
    ocp.cost.cost_type = 'LINEAR_LS'
    ocp.cost.cost_type_e = 'LINEAR_LS'

    ocp.cost.W = W
    ocp.cost.W_e = p.P
    ocp.cost.yref = np.zeros((ny,))
    ocp.cost.yref_e = np.zeros((nx,))

    # Standard stage cost
    ocp.cost.Vx = np.hstack((np.eye(nx), np.zeros((nx,nu)))).T
    ocp.cost.Vu = np.hstack((np.zeros((nu, nx)), np.eye(nu))).T
    ocp.cost.Vu[0, -1] = 1
    ocp.cost.Vx_e = np.eye(nx)

    # Set cost references (must match dimensions!)
    ocp.cost.W_0 = W.copy()
    ocp.cost.Vx_0 = np.vstack((np.eye(nx), np.zeros((nu,nx))))      # shape (5,4)
    ocp.cost.Vu_0 = np.vstack((np.zeros((nx,nu)), np.eye(nu)))      # shape (5,1)
    ocp.cost.yref_0 = np.zeros(ny)

    # Constraints
    ocp.constraints.x0 = np.zeros(4)
    ocp.constraints.lbu = np.array([p.umin])
    ocp.constraints.ubu = np.array([p.umax])
    ocp.constraints.idxbu = np.array([0])

    # Initial condition placeholder
    ocp.constraints.x0 = np.zeros(4)

    # Solver options
    ocp.solver_options.qp_solver = 'FULL_CONDENSING_QPOASES'
    ocp.solver_options.integrator_type = 'ERK'
    ocp.solver_options.nlp_solver_type = 'SQP_RTI'

    return ocp


def create_acados_ocp_cartpole2(dt, N, p):
    model = AcadosModel()
    model.name = 'cartpole_mpc'

    # === States
    x = SX.sym('x')
    v = SX.sym('v')
    th = SX.sym('th')
    w = SX.sym('w')
    x_ = vertcat(x, v, th, w)

    # === Control
    u = SX.sym('u')

    # === Parameter: u_prev (for ∆u constraint)
    u_prev = SX.sym('u_prev', 1)
    model.p = u_prev

    # === Dynamics
    f = pend_cart_nl_sym(x_, u, p)
    model.x = x_
    model.u = u
    model.f_expl_expr = f
    model.f_impl_expr = x_ - f

    # === Constraint: Δu = u - u_prev
    delta_u = u - u_prev
    model.con_h_expr = vertcat(delta_u)

    # === Setup OCP
    ocp = AcadosOcp()
    ocp.model = model
    ocp.parameter_values = np.zeros(1)

    # === Horizon
    ocp.solver_options.N_horizon = N
    ocp.solver_options.tf = dt * N

    # === Cost
    Q = np.diag([10, 10, 100, 10])
    R = np.array([[0.01]])
    W = np.block([
        [Q, np.zeros((4, 1))],
        [np.zeros((1, 4)), R]
    ])

    nx = 4
    nu = 1
    ny = nx + nu  # 4 states + 1 control
    ocp.cost.cost_type_0 = 'LINEAR_LS'
    ocp.cost.cost_type = 'LINEAR_LS'
    ocp.cost.cost_type_e = 'LINEAR_LS'

    ocp.cost.W = W
    ocp.cost.W_e = p.P  # LQR terminal cost
    ocp.cost.yref = np.zeros((ny,))
    ocp.cost.yref_e = np.zeros((nx,))

    ocp.cost.Vx = np.hstack((np.eye(nx), np.zeros((nx, nu)))).T
    ocp.cost.Vu = np.hstack((np.zeros((nu, nx)), np.eye(nu))).T
    ocp.cost.Vu[0, -1] = 1
    ocp.cost.Vx_e = np.eye(nx)

    ocp.cost.W_0 = W.copy()
    ocp.cost.Vx_0 = np.vstack((np.eye(nx), np.zeros((nu, nx))))
    ocp.cost.Vu_0 = np.vstack((np.zeros((nx, nu)), np.eye(nu)))
    ocp.cost.yref_0 = np.zeros(ny)

    # === Control input constraint
    ocp.constraints.lbu = np.array([p.umin])
    ocp.constraints.ubu = np.array([p.umax])
    ocp.constraints.idxbu = np.array([0])

    # === Δu constraint
    ocp.constraints.lh = np.array([p.dumin * dt])
    ocp.constraints.uh = np.array([p.dumax * dt])

    # === Initial state placeholder
    ocp.constraints.x0 = np.zeros(4)

    # === Solver settings
    ocp.solver_options.qp_solver = 'FULL_CONDENSING_QPOASES'
    ocp.solver_options.integrator_type = 'ERK'
    ocp.solver_options.nlp_solver_type = 'SQP_RTI'

    return ocp
