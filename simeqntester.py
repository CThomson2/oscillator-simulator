import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp, solve_bvp

n_step = 1000
dt = 0.01

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

plt.figure()
# plt.plot(positions, x[0])

# for i in range(1, n_step):
#     y0 = model(i, y0)
#     x[i] = x[i - 1] + y0[0] * dt
#     u[i] = u[i - 1] + y0[1] * dt
#     plt.plot(positions, x[i])

res = solve_ivp(model, [0, n_step], y0)


slices = []
for i in range(len(res.y[0])):
    slices.append(np.zeros(N))
    for j in range(N):
        slices[i][j] = res.y[j][i]

print(len(slices))

for i in range(len(slices)):
    if i % 5 == 0:
        plt.plot(positions, slices[i])

plt.xlabel("Position along lattice")
plt.ylabel("Oscillator displacement")
plt.show()

# res = model(0, y0)
# for i in range(N):
#     x[1][i] = x[0][i] + res[0][i] * 2
#     u[1][i] = 0 + res[1][i] * 10

# res2 = model(0, [x[1], u[1]])
# for i in range(N):
#     x[2][i] = x[1][i] + res2[0][i] * 2
#     u[2][i] = u[1][i] + res2[1][i] * 2

# res3 = model(0, [x[2], u[2]])
# for i in range(N):
#     x[3][i] = x[2][i] + res3[0][i] * 2
#     u[3][i] = u[2][i] + res3[1][i] * 2

