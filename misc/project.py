import math
import matplotlib.pyplot as plt
import time

# displacements
x = []
# no. of osc
n = 65
# tf
t = 100
# time steps (time array)
timestep = []
# string length
length = 1
# youngs mod
k = 0.1
# density
rho = 400
# alpha
a = 0
# step size inbetween osc.
h = length / (n - 1)
# mass of each osc.
m = rho * h
# stiffness
kk = k / h
# time step
dt = 0.001  # dt not h

# velocity array of time j-1
yold = []
# array to store osc. positions along string
positions = []

# print(h, m, kk)

# add values to positions array
l = 0
for i in range(n):
    l += length / n
    positions.append(l)

# velocity is zero at time t = 0
# n - 2 because boundary conditions are constant 0 at all times
for i in range(n - 2):
    yold.append(0)

# add linearly-spaced values to time array
for i in range(math.ceil(t / dt)):
    timestep.append(i)

# set initial displacement conditions for time t = 0
step = 0
for i in range(n):
    if i != 0:
        step = step + (2 * math.pi) / (n - 1)
    x.append(math.sin(step))

plt.figure(1)
plt.plot(positions, x, 'bo', positions, x, 'b')
plt.show()

# define velocity derivative function with parameters of x3 and the two balls surrounding it, x1 and x2
def models(x1, x2, x3):
    return (kk / m) * (x1 + x2 - 2 * x3) * (1 + a * (x1 - x2))

# define euler function, yold is the previous velocity value
def euler(yold, h, slope):
    return yold + h * slope

# define RK function
def RungeKutta(slope, h, x1, x2, x3, y):
    x_dum = x1
    y_dum = y
    k1 = slope

    # Second step
    x_dum = x1 + 0.5 * h
    y_dum = y + 0.5 * k1 * h
    k2 = models(x2, x3, x_dum)

    # Third step
    x_dum = x1 + 0.5 * h
    y_dum = y + 0.5 * k2 * h
    k3 = models(x2, x3, x_dum)

    # Fourth step
    x_dum = x1 + h
    y_dum = y + k3 * h
    k4 = models(x2, x3, x_dum)

    # Calculating slope
    k5 = k1 / 6 + 2 * k2 / 6 + 2 * k3 / 6 + k4 / 6

    # Calculating new y value
    return y + h * k5

# loop through time steps
for i in range(math.ceil(t / h)):
    # loop through oscillators on string
    for u in range(n - 2):
        # find slope (dv_u+1_i/dt_i+1) by calling acceleration function with current x values
        # x[u + 1] is the "current" oscillator
        slope = models(x[u + 2], x[u], x[u + 1])
        # find velocity of ball u at the next step in time using Euler with the slope that was just determined
        # u1 = euler(yold[u], h, slope)
        # same thing but with RK instead. More effective/accurate.
        u1 = RungeKutta(slope, h, x[u+1], x[u+2], x[u], yold[u])
        # set the velocity of ball u to the new velocity just found
        # note: in this program, values of v and x at "old" timesteps are simply discarded/replaced with new values
        # i.e. y and x are 1D arrays that only hold the values at the current time i
        yold[u] = u1
        # using the velocity just found (i.e. dx_u_i/dt_i) use a one-line Euler-inspired assignment to find 
        #   the ball's displacement x at that same point in time
        x[u + 1] = x[u + 1] + u1 * h

    # if (i % 5) == 0:
    #     time.sleep(2)
    #     plt.figure(2)
    #     plt.plot(le,x)
    #     plt.show

plt.plot(positions, x, 'go', positions, x, 'g')
plt.show()
# plt.savefig("Test.png",dpi=500)

# print(f"""
#       x1 = {x[2]}
#       x2 = {x[0]}
#       x3 = {x[1]}
#       kk = {kk}
#       m = {m}
#       a = {a}
#       ---------------
#       slope = {models(x[2], x[0], x[1], kk, m, a)}
#       u1    = {euler(yold[0], h, models(x[2], x[0], x[1], kk, m, a))}
#       x1    = {x[1]+euler(yold[0], h, models(x[2], x[0], x[1], kk, m, a))*h}
#       """)
