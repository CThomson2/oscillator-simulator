import math
import matplotlib.pyplot as plt
import numpy as np

L = 1
t_final = 1000
N = 65
k = 0.1
rho = 400
alpha = 0.25
pert_ini = 1.0
init_type = 1


def simulation(data):

    # ----- Parameter Definitions -----
    # x_i: displacement from equilibrium of i-th oscillator
    # u_i: velocity of i-th oscillator (eqbm vel. is zero)
    # p_i: distance along string (can be represented as i*h)
    # t: time
    # h: distance between oscillators (step size)
    # L: length of string
    # N: number of oscillators
    destruct = lambda dict, *args: (dict[arg] for arg in args)
    L, tf, N, k, rho, alpha, pert_ini, init_type = destruct(
        data, 'L', 'tf', 'N', 'k', 'rho', 'alpha', 'pert_ini', 'init_type')
    # print(L, tf, N, k, rho, alpha, pert_ini, init_type)

    h = L / (N - 1)  # step size
    m = rho * h  # mass of osc.

    def vel_model(i, t):
        return k / m * (x_rk[i] + x_rk[i - 2] - 2 * x_rk[i - 1]) * (1 + alpha * (x_rk[i] - x[i - 2]))

    u0 = 0  # x0(t) = 0
    uL = 0  # xL(t) = 0


# simulation({'L': L, 'tf': t_final, 'N': N, 'k': k, 'rho': rho, 'alpha': alpha, 'pert_ini': pert_ini, 'init_type': init_type})
