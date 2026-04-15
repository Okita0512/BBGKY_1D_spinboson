import numpy as np
from numpy import array as A


def bathParam(omega_c, alpha, ndof):
    c = A([np.sqrt(2.0) / 10])  # Model A
#   c = A([np.sqrt(2.0) / 2])   # Model B
#   c = A([np.sqrt(2.0)])       # Model C

    omega = A([1.0])

    return c, omega


class parameters:
    t = 30
    dt = 0.001
    NSteps = int(t / dt)
    NStates = 2
    M = 1
    nskip = 50

    beta, omega_c, alpha, ndof = 16, 0.0, 0.0, 1
    nmod = 1
    c, omega = bathParam(omega_c, alpha, ndof)
    ω = omega


def Hs():
    Ham = np.zeros((parameters.NStates, parameters.NStates), dtype=complex)
    Ham[0, 1] = 0.5
    Ham[1, 0] = Ham[0, 1]
    return Ham


def Qs():
    Qmat = np.zeros((parameters.ndof, parameters.NStates, parameters.NStates), dtype=complex)
    Qmat[0, 0, 0] = 1.0
    Qmat[0, 1, 1] = -1.0
    return Qmat
