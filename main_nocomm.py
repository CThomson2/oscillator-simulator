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

    destruct = lambda dict, *args: (float(dict[arg]) for arg in args)
    L, tf, N, k, rho, alpha, pert_ini = destruct(
        data, 'L', 'tf', 'N', 'k', 'rho', 'alpha', 'pert_ini')
    init_type = data['init_type']
    tf = int(tf)
    N = int(N)
    print('L', 'tf', 'N', 'k', 'rho', 'alpha', 'pert_ini', 'init_type')
    print(L, tf, N, k, rho, alpha, pert_ini, init_type)

    h = L / (N - 1)  # step size
    m = rho * h  # mass of osc.
    pos_lattice = [i * h for i in range(N)]

    dt = 1

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

    x[0] = init_disp

    # Clearly, the initial velocity is zero at every point
    init_vel = np.zeros(N)
    v[0] = init_vel

    # Initial acceleration can be found using values for displacement
    a[0][0] = 0
    a[0][N - 1] = 0
    for i in range(1, N - 1):
        a[0][i] = get_acc(0, i)

    # iterating through time-steps
    for j in range(0, n_step - 1):
        # iterating through oscillators on string
        for i in range(1, N - 1):
            v_slope = get_acc(t[j], i)
            v[j + 1][i] = v[j][i] + v_slope * dt
            x[j + 1][i] = x[j][i] + v[j][i] * dt

    # osc5 = [x[i][5] for i in range(len(t))]
    # plt.plot(t, osc5, 'r-')

    print()
    print(x[:5])
    print()
    print(v[:5])


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
