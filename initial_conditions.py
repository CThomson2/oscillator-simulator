import math
import matplotlib.pyplot as plt
import numpy as np


def get_init(L, N, amp_0, shape_0, pos):

    x_init = np.zeros(N)
    # --- Case 1 ---
    if shape_0 == 'sine':
        for i in range(N):

            x_init[i] = amp_0 * np.sin((2 * np.pi / L) * pos[i])
            # handle the potential error of boundary values being <just off> zero when they should be zero exactly
            if abs(x_init[i]) < 1e-10:
                x_init[i] = 0
    # --- Case 2 ---
    elif shape_0 == 'half-sine':
        for i in range(N):

            x_init[i] = amp_0 * np.sin((np.pi / L) * pos[i])
            if abs(x_init[i]) < 1e-10:
                x_init[i] = 0
    # --- Case 3 ---
    elif shape_0 == 'parabola':
        for i in range(N):

            x_init[i] = -(4 * amp_0 / L ** 2) * \
                (pos[i] - L / 2) ** 2 + amp_0
            if abs(x_init[i]) < 1e-10:
                x_init[i] = 0
    
    return x_init