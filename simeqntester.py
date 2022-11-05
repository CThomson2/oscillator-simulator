import math
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
from scipy.integrate import solve_ivp
from random import randint

tf = 1000
dt = 0.01
n_step = int(tf/dt)

N = 50
L = 1
h = L / (N - 1)

k = 0.1
rho = 400
alpha = 0.5

m = rho * h

init_type = 'sine'
init_amp = 1

positions = [i * h for i in range(N)]

x = []
u = []

for i in range(n_step):
    x.append( np.zeros(N) )
    u.append( np.zeros(N) )

def get_init():
    
    x_init = np.zeros(N)
    # --- Case 1 ---
    if init_type == 'sine':
        for i in range(N):

            x_init[i] = init_amp * math.sin((2 * math.pi / L) * positions[i])
            # handle the potential error of boundary values being <just off> zero when they should be zero exactly
            if abs(x_init[i]) < 1e-10:
                x_init[i] = 0
        return x_init
    # --- Case 2 ---
    elif init_type == 'half-sine':
        for i in range(N):

            x_init[i] = init_amp * math.sin((math.pi / L) * positions[i])
            if abs(x_init[i]) < 1e-10:
                x_init[i] = 0
        return x_init
    # --- Case 3 ---
    elif init_type == 'parabola':
        for i in range(N):

            x_init[i] = -(4 * init_amp / L ** 2) * \
                (positions[i] - L / 2) ** 2 + init_amp
            if abs(x_init[i]) < 1e-10:
                x_init[i] = 0
        return x_init

def get_acc(osc):
    return (k / h) / m * (osc[2] + osc[0] - 2 * osc[1]) * (1 + alpha * (osc[2] - osc[0]))

x[0] = get_init()
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

def rgb_to_hsl (r, g, b):
    r /= 255;
    g /= 255;
    b /= 255;
    l = max(r, g, b);
    s = l - min(r, g, b);

    if s:
        if l == r:
            h = (g - b) / s
        elif l == g:
            h = 2 + (b - r) / s
        else:
            h = 4 + (r - g) / s
    else:
        h = 0
        
    if 60 * h < 0:
        h = 60 * h + 360
    else:
        h = h * 60
    if s:
        if l <= 0.5:
            s = 100 * (s / (2 * l - s))
        else:
            s = 100 * (s / (2 - (2 * l - s)))
    else:
        s = 0
    l = (100 * (2 * l - s)) / 2
    
    return [h, s, l]
        
def change_col(rgb):
    p = randint(0,2)
    rgb = list(rgb)
    rgb[p] += 0.02
    if max(rgb) >= 1:
        rgb[rgb.index(max(rgb))] -= 1
    return tuple(rgb)
    

def animate(i):
    ax.clear()
    # r = (i / 100) % 1
    # g = ((100 - i) / 100) % 1
    r = abs(np.cos(0.5 * (i / 100 * 2 * np.pi + np.pi)))
    g = abs(np.sin(0.5 * (i / 100 * 2 * np.pi) + 0.2 * np.pi))
    b = abs(np.sin(0.5 * (i / 100 * 2 * np.pi - 0.75 * np.pi)))
    ax.plot(positions, slices[i], color=(r, g, b))

    ax.set_xlim([0,L])
    ax.set_ylim([-init_amp, init_amp])


ani = FuncAnimation(fig, animate, frames = n_step, interval = 25, repeat=False)


plt.xlabel("Position along lattice")
plt.ylabel("Oscillator displacement")
plt.show()
