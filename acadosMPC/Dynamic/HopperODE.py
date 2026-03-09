import numpy as np
from  scipy.integrate import odeint, solve_ivp
from Tools.Lsys import build_lsys
from Tools.mathOp import Nw, Tw
from Tools.tool import tensor_interpolation, mass_center_interpolate, ComponentForForces
from Tools.ctrl import continuous_to_discrete


class Vector(np.ndarray):

    def __new__(cls, x, u):
        # Create the ndarray instance
        obj = np.asarray(np.hstack((x, u))).view(cls)
        # Store the dimensions
        obj.nx = x.shape[0]
        obj.nu = u.shape[0]
        return obj

    def __array_finalize__(self, obj):
        # This is called whenever the array is created
        if obj is None:
            return
        self.nx = getattr(obj, 'nx', None)
        self.nu = getattr(obj, 'nu', None)

    def get_x(self):
        return self[:self.nx]

    def get_u(self):
        return self[self.nx:]

    def get(self):
        return self.get_x(), self.get_u()


def rk4(ode, t, x0, u, p=()):
    dt = t[-1] - t[0]
    k1 = ode(x0, u, p)
    k2 = ode(x0 + 0.5 * dt * k1, u, p)
    k3 = ode(x0 + 0.5 * dt * k2, u, p)
    k4 = ode(x0 + dt * k3, u, p)
    return (k1 + 2 * k2 + 2 * k3 + k4) / 6.

def solve_rk4(ode, t, x0, u, arg=()):
    y = np.zeros((x0.shape[0], t.shape[0]))
    y[:, 0] = x0
    for i, (ti, tk) in (zip(t[:-1], t[1:])):
        y[:, i+1] = y[:, i] + rk4(ode, (ti, tk), y[:,i], u, p=arg)
    return y


def solve_ode(ode_system, t, x0, args=(), method='odeint', rtol=1e-6, atol=1e-9):
    """
    Solve ODE system and return (solution, derivatives) tuple directly.

    Args:
        ode_system: HopperODE instance
        t: time points (nt,)
        x0: initial state (nx,)
        args: tuple of additional arguments (u, p)
        method: 'odeint' or 'RK45' (default: 'odeint')
        rtol: relative tolerance
        atol: absolute tolerance

    Returns:
        tuple: (solution, derivatives)
    """

    def _solve_odeint(ode, t, x0, args, rtol, atol):
        """Solve using scipy.integrate.odeint."""
        u, p = args if len(args) == 2 else (None, None)
        sol = odeint(ode, x0, t, args=(u, p), rtol=rtol, atol=atol, tfirst=True)
        ydot = np.array([ode(xi, ti, u, p) for xi, ti in zip(sol, t)])
        return sol, ydot

    def _solve_rk45(ode, t, x0, args, rtol, atol):
        """Solve using scipy.integrate.solve_ivp (RK45)."""
        u, p = args if len(args) == 2 else (None, None)
        res = solve_ivp(
            fun=lambda t, x: ode(t, x, u, p),
            t_span=(t[0], t[-1]),
            y0=x0,
            t_eval=t,
            method='RK45',
            rtol=rtol,
            atol=atol
        )
        sol = res.y.T  # Transpose to match odeint shape
        ydot = np.array([ode(ti, xi, u, p) for xi, ti in zip(sol, t)])
        return sol, ydot

    # Dispatch to the selected solver method
    solver_map = {
        'odeint': _solve_odeint,
        'RK45': _solve_rk45,
        'RK4': solve_rk4,
    }
    solver_func = solver_map.get(method, _solve_odeint)
    return solver_func(ode_system, t, x0, args, rtol, atol)


class Solver:
    """
    Enhanced solver supporting multiple integration methods from scipy.

    Supported methods:
    - 'odeint': scipy.integrate.odeint (default)
    - 'RK45': Explicit Runge-Kutta method of order 5(4)
    - 'RK23': Explicit Runge-Kutta method of order 3(2)
    - 'DOP853': Explicit Runge-Kutta method of order 8
    - 'Radau': Implicit Runge-Kutta method of Radau IIA family of order 5
    - 'BDF': Implicit multi-step variable-order (1 to 5) method
    - 'LSODA': Adams/BDF method with automatic stiffness detection
    """

    # Define which methods use solve_ivp vs odeint
    SOLVE_IVP_METHODS = {'RK45', 'RK23', 'DOP853', 'Radau', 'BDF', 'LSODA'}
    ODEINT_METHODS = {'odeint'}
    CUSTOM_METHODS = {'RK4'}

    def __init__(self, ode_system=None, t=None, x0=None, args=(), method='odeint', rtol=1e-6, atol=1e-9,
                 **solver_kwargs):
        """
        Two usage patterns:
        1. Full initialization: Solver(ode, t, x0, args) - solves immediately
        2. Partial initialization: Solver() - just creates solver object

        Args:
            ode_system: ODE system function
            t: time points
            x0: initial conditions
            args: additional arguments (u, p)
            method: integration method (see class docstring for options)
            rtol: relative tolerance
            atol: absolute tolerance
            **solver_kwargs: additional keyword arguments for solve_ivp methods
                           (e.g., max_step, first_step, dense_output, events, etc.)
        """
        self.ode = ode_system
        self.method = method
        self.rtol = rtol
        self.atol = atol
        self.solver_kwargs = solver_kwargs  # Store additional solver options
        self.x = None
        self.dx = None

        # Validate method
        if method not in (self.SOLVE_IVP_METHODS | self.ODEINT_METHODS | self.CUSTOM_METHODS):
            available_methods = sorted(self.SOLVE_IVP_METHODS | self.ODEINT_METHODS | self.CUSTOM_METHODS)
            raise ValueError(f"Unknown method '{method}'. Available methods: {available_methods}")

        # If all parameters provided, solve immediately and store results
        if ode_system is not None and t is not None and x0 is not None:
            self.x, self.dx = self._solve(t, x0, args)

    def __call__(self, ode_system, t, x0, args=(), method=None, **solver_kwargs):
        """
        Call the solver like a function.
        This RETURNS the results, doesn't store them in self.x/dx.

        Args:
            ode_system: ODE system function
            t: time points
            x0: initial conditions
            args: additional arguments (u, p)
            method: override default method for this call
            **solver_kwargs: override solver kwargs for this call
        """
        # Use provided parameters or fall back to defaults
        ode_to_use = ode_system if ode_system is not None else self.ode
        method_to_use = method if method is not None else self.method
        kwargs_to_use = {**self.solver_kwargs, **solver_kwargs}

        return self._solve_with_ode(ode_to_use, t, x0, args, method_to_use, kwargs_to_use)

    def __iter__(self):
        """
        Allows tuple unpacking: sol, deriv = solver_instance
        This only works if self.x and self.dx exist!
        """
        if self.x is not None and self.dx is not None:
            return iter((self.x, self.dx))
        else:
            raise ValueError("No stored solution. Use __call__ method or initialize with full parameters.")

    def __getitem__(self, index):
        """
        Allows indexing: sol = solver_instance[0]
        This only works if self.x and self.dx exist!
        """
        if self.x is not None and self.dx is not None:
            return (self.x, self.dx)[index]
        else:
            raise ValueError("No stored solution. Use __call__ method or initialize with full parameters.")

    def _solve(self, t, x0, args):
        """Solve using stored ode_system and method."""
        return self._solve_with_ode(self.ode, t, x0, args, self.method, self.solver_kwargs)

    def _solve_with_ode(self, ode_system, t, x0, args, method, solver_kwargs):
        """Solve with specified ode_system and method."""
        if method in self.SOLVE_IVP_METHODS:
            return self._solve_ivp(ode_system, t, x0, args, method, solver_kwargs)
        elif method in self.ODEINT_METHODS:
            return self._solve_odeint(ode_system, t, x0, args, solver_kwargs)
        elif method in self.CUSTOM_METHODS:
            return self._solve_rk4(ode_system,t, x0, args)
        else:
            raise ValueError(f"Unknown method '{method}'")

    def _solve_odeint(self, ode_system, t, x0, args, solver_kwargs):
        """Solve using scipy.integrate.odeint."""
        u, p = args if len(args) == 2 else (None, None)

        # Filter kwargs that are valid for odeint
        odeint_kwargs = {}
        valid_odeint_args = {'rtol', 'atol', 'tcrit', 'h0', 'hmax', 'hmin', 'ixpr', 'mxstep', 'mxhnil', 'mxordn',
                             'mxords'}
        for key, value in solver_kwargs.items():
            if key in valid_odeint_args:
                odeint_kwargs[key] = value

        sol = odeint(
            ode_system, x0, t,
            args=(u, p),
            rtol=self.rtol,
            atol=self.atol,
            tfirst=True,
            **odeint_kwargs
        )
        ydot = np.array([ode_system(xi, ti, u, p) for xi, ti in zip(sol, t)])
        return sol, ydot

    def _solve_ivp(self, ode_system, t, x0, args, method, solver_kwargs):
        """Solve using scipy.integrate.solve_ivp with specified method."""
        u, p = args if len(args) == 2 else (None, None)

        res = solve_ivp(
            fun=lambda t, x: ode_system(t, x, u, p),
            t_span=(t[0], t[-1]),
            y0=x0,
            t_eval=t,
            method=method,
            rtol=self.rtol,
            atol=self.atol,
            **solver_kwargs
        )

        if not res.success:
            raise RuntimeError(f"Integration failed with method '{method}': {res.message}")

        sol = res.y.T  # Transpose to match odeint shape (nt, nx)
        ydot = np.array([ode_system(ti, xi, u, p) for xi, ti in zip(sol, t)])
        return sol, ydot

    def _solve_rk4(self, ode_system, t, x0, args):
        u, p = args if len(args) == 2 else (None, None)
        x = solve_rk4(ode_system, t, x0, u, p)
        xdot = np.array([ode_system(ti, xi, u, p) for xi, ti in zip(x, t)])
        return x, xdot

    def get(self):
        """Return stored solution (backward compatibility)."""
        if self.x is not None and self.dx is not None:
            return self.x, self.dx
        else:
            raise ValueError("No stored solution. Use __call__ method or initialize with full parameters.")

    def get_available_methods(self):
        """Return list of available integration methods."""
        return sorted(self.SOLVE_IVP_METHODS | self.ODEINT_METHODS)

    def set_method(self, method):
        """Change the default integration method."""
        if method not in (self.SOLVE_IVP_METHODS | self.ODEINT_METHODS):
            available_methods = self.get_available_methods()
            raise ValueError(f"Unknown method '{method}'. Available methods: {available_methods}")
        self.method = method

    def set_tolerances(self, rtol=None, atol=None):
        """Update tolerance settings."""
        if rtol is not None:
            self.rtol = rtol
        if atol is not None:
            self.atol = atol


# USAGE EXAMPLES:

# Example 1: Different methods with immediate solving
# sol, dx = Solver(ode, tau, x0, args=(u0, param), method='RK45')
# sol, dx = Solver(ode, tau, x0, args=(u0, param), method='DOP853')
# sol, dx = Solver(ode, tau, x0, args=(u0, param), method='Radau')  # Good for stiff problems
# sol, dx = Solver(ode, tau, x0, args=(u0, param), method='BDF')    # Good for stiff problems
# sol, dx = Solver(ode, tau, x0, args=(u0, param), method='LSODA')  # Automatic stiffness detection

# Example 2: Create solver and use different methods
# solver = Solver(method='RK45', rtol=1e-8, atol=1e-10)
# sol1, dx1 = solver(ode, tau, x0, args=(u0, param))
# sol2, dx2 = solver(ode, tau, x0, args=(u0, param), method='DOP853')  # Override method for this call

# Example 3: Using additional solve_ivp options
# solver = Solver(method='RK45', max_step=0.1, dense_output=True)
# sol, dx = solver(ode, tau, x0, args=(u0, param))

# Example 4: Method-specific optimizations
# # For smooth problems - use explicit methods
# solver_smooth = Solver(method='DOP853', rtol=1e-12)  # High accuracy
#
# # For stiff problems - use implicit methods
# solver_stiff = Solver(method='Radau', rtol=1e-6, max_step=0.01)
# solver_stiff_auto = Solver(method='LSODA')  # Automatic stiffness detection
#
# # For quick exploration - use lower order methods
# solver_fast = Solver(method='RK23', rtol=1e-4)

# Example 5: Check available methods
# solver = Solver()
# print("Available methods:", solver.get_available_methods())

# Example 6: Dynamic method switching
# solver = Solver(method='RK45')
# try:
#     sol, dx = solver(ode, tau, x0, args=(u0, param))
# except RuntimeError:
#     print("RK45 failed, trying LSODA for potential stiffness...")
#     sol, dx = solver(ode, tau, x0, args=(u0, param), method='LSODA')


# DEMONSTRATION OF DIFFERENT USAGE PATTERNS:

# Pattern 1: Full initialization with immediate solving and tuple unpacking
# solver = Solver(ode, tau, x0, args=(u0, param), method="RK45")
# sol, deriv = solver  # Uses __iter__() - works because solution is stored

# Pattern 2: Full initialization with .get() method
# solver = Solver(ode, tau, x0, args=(u0, param), method="RK45")
# sol, deriv = solver.get()  # Uses stored solution

# Pattern 3: Partial initialization, then __call__
# solver = Solver(method="RK45")
# sol, deriv = solver(ode, tau, x0, args=(u0, param))  # Uses __call__() - returns directly

# Pattern 4: Create empty solver, use as callable
# solver = Solver()
# sol, deriv = solver(ode, tau, x0, args=(u0, param), method="RK45")  # Uses __call__()


# SIMPLIFIED VERSION - What you probably want:
class SimpleSolver:
    """
    Simplified solver that supports both initialization patterns.
    """

    def __init__(self, ode_system=None, t=None, x0=None, args=(), method='odeint', rtol=1e-6, atol=1e-9):
        self.method = method
        self.rtol = rtol
        self.atol = atol
        self.x = None
        self.dx = None

        # If all required params provided, solve immediately
        if ode_system is not None and t is not None and x0 is not None:
            self.x, self.dx = self._solve_ode(ode_system, t, x0, args)

    def __iter__(self):
        """Enable: sol, deriv = SimpleSolver(ode, t, x0, args)"""
        if self.x is None:
            raise ValueError("No solution stored. Provide all parameters in __init__.")
        return iter((self.x, self.dx))

    def __call__(self, ode_system, t, x0, args=()):
        """Enable: sol, deriv = solver(ode, t, x0, args)"""
        return self._solve_ode(ode_system, t, x0, args)

    def _solve_ode(self, ode_system, t, x0, args):
        """Core solving logic."""
        u, p = args if len(args) == 2 else (None, None)

        if self.method == 'RK45':
            res = solve_ivp(
                fun=lambda t, x: ode_system(t, x, u, p),
                t_span=(t[0], t[-1]),
                y0=x0,
                t_eval=t,
                method='RK45',
                rtol=self.rtol,
                atol=self.atol
            )
            sol = res.y.T
        else:  # odeint
            sol = odeint(ode_system, x0, t, args=(u, p), rtol=self.rtol, atol=self.atol, tfirst=True)

        ydot = np.array([ode_system(ti, xi, u, p) for xi, ti in zip(sol, t)])
        return sol, ydot


# USAGE EXAMPLES:
# Method 1 - Direct unpacking (uses __iter__):
# sol, deriv = SimpleSolver(ode, tau, x0, args=(u0, param), method="RK45")

# Method 2 - Callable style (uses __call__):
# solver = SimpleSolver(method="RK45")
# sol, deriv = solver(ode, tau, x0, args=(u0, param))

# Usage examples:
# Option A - Direct function call:
# sol, deriv = solve_ode(ode, tau, x0, args=(u0, param), method="RK45")

# Option B - Class with tuple unpacking:
# sol, deriv = Solver(ode, tau, x0, args=(u0, param), method="RK45")

# Option C - Original style (still works):
# solver = Solver(ode, tau, x0, args=(u0, param), method="RK45")
# sol, deriv = solver.get()

def jacobian_HopperODE(x, u, p):
    return build_lsys(x, u, p)


def jacobian_HopperODE_discrete(x, u, p):
    Ac, Bc = jacobian_HopperODE(x, u, p)
    Ad, Bd = continuous_to_discrete(Ac, Bc, p.sol.dt)
    return Ad, Bd


class HopperODE:

    def __init__(self, nx, nu):
        self.nx = nx
        self.nu = nu

    def _ode(self, x, u, p):
        """
        NL_HOPPERDYN Summary of this function goes here
        Detailed explanation goes here

        Parameters:
        x(13) - State Vector 6DoF: px,py,pz, vx,vy,vz, roll, pitch,yaw,
              - wx, wy, wz, m
        u(9)  - Control Vector: Fx,Fy,Fz, Mx,My,Mz
        param - System Parameters: Jxx,Jyy,Jzz,Jxy,Jyz,Jxz,m,Cx,Cy,Cz,Cr,Cq,Cp,g,
              - beta, gamma, rho, h, d, k_quat
        """

        # State vector
        # r = x[0:3]           # position
        v = x[3:6]  # velocity
        a = x[6:9]  # Euler's angles
        w = x[9:12]  # angular velocity
        m = x[12]  # current mass

        # Disturbances
        Fd = p.d[0:3]  # Disturbance Forces
        Md = p.d[3:6]  # Disturbance Moments

        # Mass and Inertia matrices
        mdot = -np.abs(u[0]) * p.beta
        rcg = mass_center_interpolate(x[12], p.mass, p.rcgs)
        Hf, Hm = ComponentForForces(p.geom, rcg)

        # Control Forces
        uuu = np.vstack([Hf, Hm]) @ u
        F = uuu[0:3]  # Forces
        M = uuu[3:6]  # Moments
        Jo = tensor_interpolation(m, p.mass, p.Inertia)
        Jdot = p.Jdot_mlt * mdot  # debug

        Mass = np.diag([m, m, m])

        # Damping matrices
        Cv = p.Cv
        Cw = p.Cw

        # Gravity
        G = p.gc * Nw(a) @ p.G    # Get gravity in body frame

        # System of Equations
        xdot = np.hstack(
            (
                v,
                -np.cross(w, v) + np.linalg.solve(Mass, (-Cv @ v + F + Fd)) + G,
                Tw(a) @ w,
                np.linalg.solve(Jo, (-np.cross(w, Jo @ w) - Cw @ w - Jdot @ w + M + Md)),
                mdot
            )
        )

        return xdot


    def __call__(self, t, x, *args):
        u = args[0]
        p = args[1]
        return self._ode(x, u, p)


class HopperActDyn(HopperODE):

    def _actuator_dynamic(self, u, u_app, p):
        # Acting dynamic is T * u_dot + t1 * u = u_smc
        # T = p.T
        # t1 = p.t1
        # u - current or previous control magnitude
        # u_smc - applied control
        du_dt = (u_app - p.t1 * u) / p.T
        return du_dt

    def _ode(self, s, u, p):
        x_ = s[:self.nx]
        u_ = s[self.nx:]
        return np.hstack((super()._ode(x_, u_, p),
                          self._actuator_dynamic(u_, u, p.ctrl)
                          )
                         )
