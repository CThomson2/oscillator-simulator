import math
import matplotlib.pyplot as plt
import numpy as np


def get_init(L, N, pert_ini, init_type, h, pos_lattice):
    # Define arrays to store initial velocity at every oscillator i
    v_init = np.zeros(N)
    # --- Case 1 ---
    if init_type == 'sine':
        for i in range(N):

            v_init[i] = pert_ini * math.sin((2 * math.pi / L) * pos_lattice[i])
            # handle the potential error of boundary values being <just off> zero when they should be zero exactly
            if abs(v_init[i]) < 1e-10:
                v_init[i] = 0
        return v_init
    # --- Case 2 ---
    elif init_type == 'half-sine':
        for i in range(N):

            v_init[i] = pert_ini * math.sin((math.pi / L) * pos_lattice[i])
            if abs(v_init[i]) < 1e-10:
                v_init[i] = 0
        # plt.plot(lattice, v_init, 'r.-')
        # plt.show()
        return v_init
    # --- Case 3 ---
    elif init_type == 'parabola':
        for i in range(N):

            v_init[i] = -(4 * pert_ini / L ** 2) * \
                (i * h - L / 2) ** 2 + pert_ini
            if abs(v_init[i]) < 1e-10:
                v_init[i] = 0
        # plt.plot(lattice, v_init, 'b.-')
        # plt.show()
        return v_init
    # Add more polynomial possibilites here (cubic, quartic etc.)
