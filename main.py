# import packages
import math
import matplotlib.pyplot as plt
import numpy as np
import webbrowser
import os

# import modules
from initial_conditions import get_init

# Boundary Conditions
# since x_0(t) = x_L(t) = 0, velocity at the boundaries is also a constant zero
v0 = 0
vL = 0

# Move this function to its own file later


# def euler(t, x):


def simulate(data):

    # ----- Parameter Definitions -----
    # x_i: displacement from equilibrium of i-th oscillator
    # u_i: velocity of i-th oscillator (eqbm vel. is zero)
    # p_i: distance along string (can be represented as i*h)
    # t: time
    # h: distance between oscillators (step size)
    # L: length of string
    # N: number of oscillators
    destruct = lambda dict, *args: (float(dict[arg]) for arg in args)
    L, tf, N, k, rho, alpha, pert_ini = destruct(
        data, 'L', 'tf', 'N', 'k', 'rho', 'alpha', 'pert_ini')
    init_type = data['init_type']
    N = int(N)
    print(L, tf, N, k, rho, alpha, pert_ini, init_type)

    # calculate stepsize (distance between oscillators) and mass of each oscillator
    h = L / (N - 1)  # step size
    m = rho * h  # mass of osc.
    # array holding positions of oscillators along string
    pos_lattice = [i * h for i in range(N)]
    # choose a timestep
    dt = 0.01
    n_step = math.ceil(tf / dt)

    # def vel_model(i, t):
    #     return k / m * (x_rk[i] + x_rk[i - 2] - 2 * x_rk[i - 1]) * (1 + alpha * (x_rk[i] - x[i - 2]))

    # -----
    # Aim: to find the velocity of each oscillator at each point in time
    # Method: create a 2D array to store arrays of velocities at each time-step for each oscillator
    # -----

    # STEP 1: create function to find initial displacement from user-specified conditions
    # > function exported to and located in initial_conditions.py
    init_disp = get_init(L, N, pert_ini, init_type, h, pos_lattice)
    plt.plot([i * h for i in range(N)], init_disp, 'r.-')
    plt.savefig(f'./static/img/graph.png')
    # plt.show()

    # -----
    # show user the initial conditions in graphical form
    # filename = 'file:///' + os.getcwd() + '/' + 'initcond.html'
    # webbrowser.open_new_tab(filename)
    # -----

    # -----
    # Clearly, the initial velocity is zero at every point
    init_vel = np.zeros(N)
    # -----

    # STEP 2:

    # for i in range(N):

    # t_eul_i = np.zeros(n_step + 1)
    # v_eul_i = np.zeros(n_step + 1)
    # t_eul_i[0] = 0
    # v_eul_i[0] =
    # for j in range(n_step):
    #     v_eul_i[j + 1] = euler(j, i)

# -----


# Hardcoded parameters for testing and debugging
L = 1
t_final = 10
N = 8
k = 0.1
rho = 400
alpha = 0.25
pert_ini = 3.0
init_type = 'sine'

simulate({'L': L, 'tf': t_final, 'N': N, 'k': k, 'rho': rho,
         'alpha': alpha, 'pert_ini': pert_ini, 'init_type': init_type})
