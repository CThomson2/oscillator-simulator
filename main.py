# import packages
import math
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
from scipy.integrate import solve_ivp
from random import randint
from initial_conditions import get_init

def Fermi(data):

    destruct_dict = lambda dict, *args: (float(dict[arg]) for arg in args)

    L, tf, N, k, rho, alpha, amp_0 = destruct_dict(data, 'L', 'tf', 'N', 'k', 'rho', 'alpha', 'amp_0')
    shape_0 = data['shape_0']

    N = int(N)
    tf = int(tf)

    dt = 0.01
    n_step = int(tf/dt)
    h = L / (N - 1)
    m = rho * h
    positions = [i * h for i in range(N)]

    x = []
    u = []

    for i in range(n_step):
        x.append( np.zeros(N) )
        u.append( np.zeros(N) )

    def get_acc(osc):
        return (k / h) / m * (osc[2] + osc[0] - 2 * osc[1]) * (1 + alpha * (osc[2] - osc[0]))

    x[0] = get_init(L, N, amp_0, shape_0, positions)
    u[0] = np.zeros(N)

    def model(t, y):
        x = y[:int(len(y) / 2)]
        u = y[int(len(y) / 2):]

        fx = u
        fu = np.zeros(N)

        for i in range(1, N-1):
            fu[i] = get_acc([ x[i - 1], x[i], x[i + 1] ])

        return np.concatenate([fx, fu])

    y0 = np.concatenate([x[0], u[0]])

    t_eval = np.linspace(0, tf, num=n_step)
    res = solve_ivp(model, [0, tf], y0)

    slices = []
    for i in range(len(res.y[0])):
        slices.append(np.zeros(N))
        for j in range(N):
            slices[i][j] = res.y[j][i]

    fig, ax = plt.subplots()

    print(len(slices))
           
    def animate(i):
        ax.clear()
        
        # The following lines change the graph colour over time
        # Each RGB channel follows a different sinusoid, with the peaks of all
        # three seperated
        # The amplitude is in order to prevent unaesthethic light colours
        r = abs(0.75*np.cos(0.5 * (i / 100 * 2 * np.pi + np.pi)))     
        g = abs(0.75*np.sin(0.5 * (i / 100 * 2 * np.pi) + 0.2 * np.pi))
        # b = abs(0.75*np.sin(0.5 * (i / 100 * 2 * np.pi - 0.75 * np.pi)))
        b = 0.7
        ax.plot(positions, slices[i], marker='o', markersize=3, color=(r, g, b))

        ax.set_xlim([0,L])
        ax.set_ylim([-amp_0, amp_0])


    ani = FuncAnimation(fig, animate, frames = n_step - 1, interval = 1, repeat=False)

    plt.title("Oscillator displacement along lattice (string) over time")
    plt.xlabel("Position along lattice")
    plt.ylabel("Oscillator displacement")
    plt.show()


tf = 2500

N = 100
L = 1

k = 0.4
rho = 200
alpha = 0.9

amp_0 = 1
shape_0 = 'sine'

Fermi({'L': L, 'tf': tf, 'N': N, 'k': k, 'rho': rho, 'alpha': alpha, 'amp_0': amp_0, 'shape_0': shape_0})