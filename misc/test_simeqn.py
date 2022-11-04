import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.integrate import solve_ivp

L = 1
N = 6
k = 0.1
rho = 400
h = L / (N - 1)
m = rho * h
alpha = 0
pert = 1
pos = [i * h for i in range(N)]

def get_init(N = N, pert = pert, L = L):
    xin = np.zeros(N)
    for i in range(N):
        xin[i] = pert * math.sin((2 * math.pi / L) * pos[i])
        if abs(xin[i]) < 1e-10:
            xin[i] = 0
    return xin

def get_acc(osc):
    result = (k / h) / m * (osc[2] + osc[0] - 2 * osc[1]) * (1 + alpha * (osc[2] - osc[0]))
    return result

def model(t, y):
    x0 = 0
    v0 = 0
    x5 = 0
    v5 = 0
    x1 = y[0]
    x2 = y[1]
    x3 = y[2]
    x4 = y[3]
    v1 = y[4]
    v2 = y[5]
    v3 = y[6]
    v4 = y[7]
    f1 = v1
    f2 = v2
    f3 = v3
    f4 = v4
    f5 = get_acc(
        [x0, x1, x2]
    )
    f6 = get_acc(
        [x1, x1, x3]
    )
    f7 = get_acc(
        [x2, x3, x4]
    )
    f8 = get_acc(
        [x3, x4, x5]
    )
    return [f1, f2, f3, f4, f5, f6, f7, f8]

# plt.plot(get_init())
# plt.show()

xin = get_init()
initial = np.concatenate([xin[1:-1], np.zeros(4)])
tfinal = 30

t_eval = np.linspace(0, tfinal, num=500)
x = solve_ivp(model, [0, tfinal], initial, t_eval=t_eval)

print(len(x.y))

plt.plot(range(6), [0, x.y[0][10], x.y[1][10], x.y[2][10], x.y[3][10], 0])
plt.plot(range(6), [0, x.y[0][20], x.y[1][20], x.y[2][20], x.y[3][20], 0])
plt.plot(range(6), [0, x.y[0][30], x.y[1][30], x.y[2][30], x.y[3][30], 0])
plt.plot(range(6), [0, x.y[0][40], x.y[1][40], x.y[2][40], x.y[3][40], 0])
plt.plot(range(6), [0, x.y[0][50], x.y[1][50], x.y[2][50], x.y[3][50], 0])
plt.plot(range(6), [0, x.y[0][60], x.y[1][60], x.y[2][60], x.y[3][60], 0])
plt.plot(range(6), [0, x.y[0][70], x.y[1][70], x.y[2][70], x.y[3][70], 0])
plt.plot(range(6), [0, x.y[0][80], x.y[1][80], x.y[2][80], x.y[3][80], 0])
plt.plot(range(6), [0, x.y[0][90], x.y[1][90], x.y[2][90], x.y[3][90], 0])
plt.plot(range(6), [0, x.y[0][100], x.y[1][100], x.y[2][100], x.y[3][100], 0])
plt.plot(range(6), [0, x.y[0][110], x.y[1][110], x.y[2][110], x.y[3][110], 0])
plt.plot(range(6), [0, x.y[0][120], x.y[1][120], x.y[2][120], x.y[3][120], 0])
plt.plot(range(6), [0, x.y[0][130], x.y[1][130], x.y[2][130], x.y[3][130], 0])
plt.plot(range(6), [0, x.y[0][140], x.y[1][140], x.y[2][140], x.y[3][140], 0])
plt.plot(range(6), [0, x.y[0][150], x.y[1][150], x.y[2][150], x.y[3][150], 0])
plt.plot(range(6), [0, x.y[0][160], x.y[1][160], x.y[2][160], x.y[3][160], 0])
plt.plot(range(6), [0, x.y[0][170], x.y[1][170], x.y[2][170], x.y[3][170], 0])
plt.plot(range(6), [0, x.y[0][180], x.y[1][180], x.y[2][180], x.y[3][180], 0])
plt.plot(range(6), [0, x.y[0][190], x.y[1][190], x.y[2][190], x.y[3][190], 0])
plt.plot(range(6), [0, x.y[0][200], x.y[1][200], x.y[2][190], x.y[3][190], 0])

plt.show()