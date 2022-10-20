# import packages
import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
import webbrowser
import os

# import modules
from initial_conditions import get_init

# Boundary Conditions
# since x_0(t) = x_L(t) = 0, velocity at the boundaries is also a constant zero
x0 = 0
xL = 0
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
    tf = int(tf)
    N = int(N)
    print('L', 'tf', 'N', 'k', 'rho', 'alpha', 'pert_ini', 'init_type')
    print(L, tf, N, k, rho, alpha, pert_ini, init_type)

    # calculate stepsize (distance between oscillators) and mass of each oscillator
    h = L / (N - 1)  # step size
    m = rho * h  # mass of osc.
    # array holding positions of oscillators along string
    # will be used as the x-axis for plotting variables of interest along length of string
    pos_lattice = [i * h for i in range(N)]
    # choose a timestep (also in ms)

    # IMPORTANT - dt shouldn't be hardset to 1. Sometimes that will be too large.
    dt = 1
    # But, if it's less than 1, the t array will contain decimals (e.g. t = [0, 0.2, 0.4, ... , 99.8, 100])
    # Clearly, you'll need to normalise the times so that the time step is equal to 1, which will increase
    # the final time. This means that when you plot the time, you need to use the original timescale provided
    # by the user.

    n_step = math.ceil(tf / dt)
    # create time-step array
    t = np.linspace(0, tf - 1, n_step, dtype=int)
    # create arrays that will hold all displacements, velocities and accelerations in the 2D space of discretised time and distance along string
    x = []
    v = []
    a = []
    # 2D array will hold an array of displacements for each time-step
    for i in range(n_step):
        x.append(np.zeros(N))
        v.append(np.zeros(N))
        a.append(np.zeros(N))

    # define governing acceleration equation for ease of access
    def get_acc(t, i):
        # check for boundary conditions
        if i in [0, N - 1]:
            return 0

        # otherwise use i parameter to determine the acceleration of the i-th oscillator
        return k / m * (x[t][i + 1] + x[t][i - 1] - 2 * x[t][i]) * (1 + alpha * (x[t][i + 1] - x[t][i - 1]))

    # def vel_model(i, t):
    #     return k / m * (x_rk[i] + x_rk[i - 2] - 2 * x_rk[i - 1]) * (1 + alpha * (x_rk[i] - x[i - 2]))

    # -----
    # Aim: to find the velocity of each oscillator at each point in time
    # Method: create a 2D array to store arrays of velocities at each time-step for each oscillator
    # -----

    # -----
    # STEP 1: create function to find initial displacement from user-specified conditions
    # Additionally, store initial velocities and accelerations in main solution array
    # -----

    # function get_init exported to and located in initial_conditions.py
    init_disp = get_init(L, N, pert_ini, init_type, h, pos_lattice)
    # plt.plot([i * h for i in range(N)], init_disp, 'r.-')
    # plt.savefig(f'./static/img/graph.png')
    # plt.show()

    # -----
    # show user the initial conditions in graphical form
    # filename = 'file:///' + os.getcwd() + '/' + 'initcond.html'
    # webbrowser.open_new_tab(filename)
    # -----

    #

    x[0] = init_disp

    # Clearly, the initial velocity is zero at every point
    init_vel = np.zeros(N)
    v[0] = init_vel

    # Initial acceleration can be found using values for displacement
    a[0][0] = 0
    a[0][N - 1] = 0
    for i in range(1, N - 1):
        a[0][i] = get_acc(0, i)

    # plt.figure(figsize=(6, 3))
    # plt.subplot(131)
    # plt.plot(pos_lattice, x[0], 'b.-')
    # plt.xlabel('string')
    # plt.ylabel('displacement')
    # plt.subplot(132)
    # plt.plot(pos_lattice, a[0], 'r.-')
    # plt.xlabel('string')
    # plt.ylabel('acceleration')

    # plt.show()

    # print('Mass:', m, 'Step size h:', h, 'Length density rho:', rho)
    # print(
    #     f'Init. Disp: {x[0]}, \n\nInit. Vel: {v[0]}, \n\nInit. Acc: {a[0]}\n')

    # -----
    # STEP 2: Taking slices of time, we study the velocities and displacements of each oscillator using provided recurrence relations
    # -----

    # -----
    # To find displacement at any given point, we need its derivative: velocity.
    # To find velocity, we need the acceleration at that point.
    # Thus, we need to fill the velocity and displacement arrays iteratively.
    # i.e. we find velocity at point i at time-step j using the formula for acceleration at point i at time-step j-1
    # We do this for all oscillators i in-between the boundary conditions, which remain at zero displacement throughout
    # Next, we find displacement at point i at time-step j using the indexed value for velocity that we just found and
    #   stored in the appropriate location in the global velocity array.
    # Before this algorithm, we can fill some rows in the 2D grids (arrays) for v and x
    # We have already inserted the initial conditions, but we also know the displacement at time 0 + ∂t (i.e. first time step)
    # Since v_0(t) = 0, the displacement cannot change in the second time-step. Therefore x_i(t1) = x_i(t0)

    # print(f'\n{t}\n')

    # throughout this program, i is associated with position along string and j with time-step
    # pattern is continued here for clarity and consistency
    # iterating through time-steps
    for j in range(0, n_step - 1):
        # iterating through oscillators on string
        for i in range(1, N - 1):
            v_slope = get_acc(t[j], i)
            v[j + 1][i] = v[j][i] + v_slope * dt
            x[j + 1][i] = x[j][i] + v[j][i] * dt

    # b1 = {
    #     "x": [],
    #     "v": []
    # }
    # for ball in x:
    #     b1["x"].append(ball[4])
    # for ball in v:
    #     b1["v"].append(ball[4])

    # plt.plot(t, b1["x"], 'r.-', t, b1["v"], 'g.-')

    # plt.figure(figsize=(9, 6))

    # plt.subplot(131)
    # plt.plot(pos_lattice, x[0])
    # plt.subplot(132)
    # plt.axes((0, -1.0, 1.0, 2.0))
    # plt.plot(pos_lattice, x[100])
    # # plt.axes((0, -1.0, 1.0, 1.0))
    # plt.subplot(133)
    # plt.plot(pos_lattice, x[200])
    # plt.axes((0, -1.0, 1.0, 1.0))
    # plt.subplot(1, 3, (2, 1))
    # plt.plot(pos_lattice, x[3])
    # plt.subplot(1, 3, (2, 2))
    # plt.plot(pos_lattice, x[4])

    plt.plot(pos_lattice, x[0], 'b-')
    # plt.plot(pos_lattice, x[50], 'r--')
    plt.plot(pos_lattice, x[80], 'b.-')
    plt.plot(pos_lattice, x[1180], 'g-')
    # plt.plot(pos_lattice, x[2000], 'g--')
    # plt.plot(pos_lattice, x[2500], 'g.-')
    # plt.plot(pos_lattice, x[3000], 'r-')
    # plt.plot(pos_lattice, x[3500], 'r--')
    # plt.plot(pos_lattice, x[7000], 'r.-')

    # osc5 = [x[i][5] for i in range(len(t))]
    # plt.plot(t, osc5, 'r-')

    plt.show()

    # Hardcoded parameters for testing and debugging


L = 1
t_final = 5000
N = 60
k = 0.1
rho = 400
alpha = 0
pert_ini = 1.0
init_type = 'parabola'

simulate({'L': L, 'tf': t_final, 'N': N, 'k': k, 'rho': rho,
         'alpha': alpha, 'pert_ini': pert_ini, 'init_type': init_type})
