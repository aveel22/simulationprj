import numpy as np
from scipy.interpolate import interp1d
from Tools.mathOp import Nw, N2B, skew, rotx, roty, rotz


def JetForceOnComponents(parameters, rc):
    # JETFORCEONCOMPONENTS Summary of this function goes here
    # Detailed explanation goes here

    eul_f = -np.deg2rad(parameters[3:6])
    eul_L = -np.deg2rad(parameters[6:9])
    hf = Nw(eul_f) @ np.array([1, 0, 0])
    L = Nw(eul_L) @ np.asarray([0,  parameters[1], parameters[2]])
    r = L + np.array([parameters[0], 0, 0])
    hm = skew(r - rc.T) @ hf
    return hf, hm


def ForceOnComponents(parameters, rc):
    # %FORCEONCOMPONENTS Summary of this function goes here
    # %   Detailed explanation goes here
    # % parameters - geomentry
    # % rc = [xc, yc, zc] - center of mass

    # % Single force

    eul_f = -np.deg2rad(parameters[3:6])
    eul_L = -np.deg2rad(parameters[6:9])
    hf = Nw(eul_f) @ np.asarray([1,0,0])
    L = Nw(eul_L) @ np.asarray([parameters[0], 0, 0])
    r = L + np.asarray([parameters[1], 0, 0])
    hm = skew(r - rc.T) @ hf
    return hf, hm


def  ComponentForForces(geom, rc, skip_column=None):
    # %COMPONENTFORFORCES Summary of this function goes here
    # %   Detailed explanation goes here

    if skip_column is None:
        skip_column = []

    Hf = np.zeros((3,9))
    Hm = np.zeros((3,9))
    k = 0
    for i in range(9):
        if i in skip_column:
            continue

        if i == 0:
            [hf, hm] = JetForceOnComponents(geom[i,:], rc)
        else:
            [hf, hm] = ForceOnComponents(geom[i,:], rc)

        Hf[:, k] = hf
        Hm[:, k] = hm
        k += 1
    return Hf, Hm

def  ComponentForForces_ENU(geom, rc, C2B=None, skip_column=None):
    # %COMPONENTFORFORCES Summary of this function goes here
    # %   Detailed explanation goes here
    # geom:     array of coordinate of actuator forces
    # rc:       array of centers of masses
    # eul:      angle rotation
    #           eul = np.array([0, -0.5 * np.pi, 0.5 * np.pi]) if N2B
    #           eul = np.array([0, 0.5 * np.pi, 0]) if Nw

    if skip_column is None:
        skip_column = []
    if C2B is None:
        C2B = np.eye(len(rc))

    xc = C2B.T @ rc
    Hf = np.zeros((rc.shape[0],geom.shape[0]))
    Hm = np.zeros_like(Hf)
    k = 0
    for i in range(geom.shape[0]):
        if i in skip_column:
            continue

        if i == 0:
            [hf, hm] = JetForceOnComponents(geom[i,:], xc)

        else:
            [hf, hm] = ForceOnComponents(geom[i,:], xc)

        Hf[:, k] = C2B @ hf
        Hm[:, k] = C2B @ hm
        k += 1
    return Hf, Hm


def Vector2Tensor(vector):
    return np.array([
        [vector[0], vector[3], vector[4]],
        [vector[3], vector[1], vector[5]],
        [vector[4], vector[5], vector[2]],
    ])

def Tensor2Vector(tensor):
    return np.array([tensor[0,0], tensor[1,1], tensor[2,2], tensor[0,1], tensor[0,2], tensor[1,2]])


def  RotateTensor(J, R=None, eul=None):
    # J - tensor of inertia
    # R - transition matrix
    # eul - Euler's angles: to support legacy

    if eul is not None:
        # UnpackEuler's angles
        phi = eul[0]
        theta = eul[1]
        psi = eul[2]
        # %%Build Rotation Matrix
        R = rotz(theta) @ roty(psi) @ rotx(phi)

    # Tensor Rotation
    J_rotated = R @ J @ R.T
    return J_rotated

def DiagTensorInertia(tensor):
    return np.diag(np.diag(tensor))


def tensor_interpolation(t_input, mass, Tensors):
    """
    Interpolates a 3D tensor along its third dimension using spline interpolation.

    Parameters:
        t_input : float
            Input mass value where interpolation is desired
        mass : array_like
            1D array of mass values corresponding to tensor slices
        Tensors : ndarray
            3D array of shape (N, M, K) where K matches length(mass)

    Returns:
        T_interp : ndarray
            Interpolated 2D tensor of shape (N, M)
    """
    N, M, K = Tensors.shape

    # Preallocate the output tensor
    T_interp = np.zeros((N, M))
    kind = 'cubic'if K > 2 else 'linear'

    # Perform interpolation for each element
    for i in range(N):
        for j in range(M):
            # Create interpolation function for this component
            f_interp = interp1d(mass, Tensors[i, j, :], kind=kind)
            # Evaluate at desired point
            T_interp[i, j] = f_interp(t_input)

    return T_interp


def tensor_interpolation_vectorized(t_input, mass, Tensors):
    # Advanced alternative
    # Reshape tensor to 2D (N*M, K) and interpolate all at once
    reshaped = Tensors.reshape(-1, len(mass))
    interp_values = np.array([interp1d(mass, row, kind='cubic')(t_input)
                          for row in reshaped])
    return interp_values.reshape(Tensors.shape[0], Tensors.shape[1])


def mass_center_interpolate(m_input, mass, rcgs):
    """
    Interpolates the center of mass coordinates based on current mass.

    Parameters:
        m_input : float
            Current mass value where interpolation is desired
        mass : array_like
            1D array of mass values
        rcgs : ndarray
            2D array of center of mass coordinates (N x 3) where N matches len(mass)

    Returns:
        rcg_interp : ndarray
            Interpolated center of mass coordinates [x, y, z]
    """
    # Initialize output array
    rcg_interp = np.zeros(3)

    # Interpolate each coordinate component
    for j in range(3):
        # Create interpolation function for this coordinate
        f_interp = interp1d(mass, rcgs[:, j], kind='cubic')
        # Evaluate at desired mass
        rcg_interp[j] = f_interp(m_input)

    return rcg_interp


def mass_center_interpolate_vectorized(m_input, mass, rcgs):
    """Vectorized version using interpolation for all axes at once"""
    interp_funcs = [interp1d(mass, rcgs[:, j], kind='cubic') for j in range(3)]
    return np.array([f(m_input) for f in interp_funcs])