# import packages
import math
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
from scipy.integrate import solve_ivp
from scipy.fft import fft
from random import randint
from initial_conditions import get_init
import time

def Fermi(data):

    # create a function that breaks the dictionary of user inputs down into its constituent parameters
    destruct_dict = lambda dict, *args: (float(dict[arg]) for arg in args)

    # extract all variables from the data dictionary
    L, tf, N, k, rho, alpha, amp_0 = destruct_dict(data, 'L', 'tf', 'N', 'k', 'rho', 'alpha', 'amp_0')
    shape_0 = data['shape_0']

    # set following params to integer format so they can be iterated through later
    N = int(N)
    tf = int(tf)

    dt = 1
    # calculate total time steps the variables of interest will be recorded at
    n_step = int(tf/dt)
    # calculate stepsize (distance between oscillators) and mass of each oscillator
    h = L / (N - 1)
    m = rho * h
    # create array holding positions of oscillators along string
    # will be used as the x-axis for plotting variables of interest along length of string
    positions = [i * h for i in range(N)]

    # create arrays that will hold all displacements and velocities in the
    #   2D space of discretised time and position on the string
    x = []
    u = []

    # add empty array to hold values for each slice in time
    for i in range(n_step):
        x.append( np.zeros(N) )
        u.append( np.zeros(N) )

    # define governing acceleration formula with a list parameter of relevant oscillator displacements
    def get_acc(osc):
        return (k / h) / m * (osc[2] + osc[0] - 2 * osc[1]) * (1 + alpha * (osc[2] - osc[0]))

    # set the displacements and velocities at t = 0 to the initial conditions derived from user input
    #   and calculated in seperate file
    x[0] = get_init(L, N, amp_0, shape_0, positions)
    u[0] = np.zeros(N)

    # define model function to be used in solve_ivp later
    # this function returns the derivatives (slopes) of all input variables
    def model(t, y):
        # extract the displacements and velocities from y array
        x = y[:int(len(y) / 2)]
        u = y[int(len(y) / 2):]

        # the derivative of displacement is simply the velocity
        fx = u
        fu = np.zeros(N)

        # fill dv/dt array using acceleration function and provided displacement values
        for i in range(1, N-1):
            fu[i] = get_acc([ x[i - 1], x[i], x[i + 1] ])

        return np.concatenate([fx, fu])

    # set the initial model input to the combined initial conditions
    y0 = np.concatenate([x[0], u[0]])

    # use scipy's solve_ivp method to solve the Fermi-Pasta problem using Runge-Kutta (4th order)
    t_eval = np.linspace(0, tf, num=n_step)
    res = solve_ivp(model, [0, tf], y0, t_eval=t_eval)

    # solve_ivp returns an multidimensional array of arrays, each nested array holding all
    #   displacements for a specific oscillator
    # to graph the string through time, the provided array must be transposed so that it will
    #   hold nested arrays of time-dependent strings (where each string represenets the displacements
    #   of its oscillators).
    slices = []
    for i in range(len(res.y[0])):
        slices.append(np.zeros(N))
        for j in range(N):
            slices[i][j] = res.y[j][i]


    # FOURIER ANALYSIS

    # create array of fourier coefficients for the N oscillators at every time step
    fourier = []
    for i in range(0, len(slices)):
        fourier.append(np.absolute(np.fft.rfft(slices[i])) / (N / 2))

    # normalise fourier coefficients at any given time step
    def normalize(vals):
        total = sum(vals)
        return [v / total for v in vals]

    # for i in range(len(fourier)):
    #     fourier[i] = normalize(fourier[i])

    frequencies = []
    for j in range(len(fourier[1])):
        frequencies.append(np.zeros(len(fourier)))
        for i in range(len(fourier)):
            frequencies[j][i] = fourier[i][j]


    plt.figure(figsize=(10, 5))
    plt.plot(frequencies[1], color='blue', label='First')
    plt.plot(frequencies[2], color='orange', label='Second')
    plt.plot(frequencies[3], color='green', label='Third')
    plt.plot(frequencies[4], color='red', label='Fourth')
    plt.legend()

    plt.xlim([0, tf])
    plt.xlabel("time, t [ms]")
    plt.ylabel("Fourier Coefficients")
    plt.title("Fourier series of lattice oscillations at discrete time steps")

    plt.show()

    # the final stage is plotting the results
    # using MatPlotLib's animation class, the solution is plotted on a live graph

    # fig, ax = plt.subplots()
    # fig.set_size_inches(12, 6)

    # def animate(i):
    #     # reset the graph at each step
    #     ax.clear()
        
    #     # the following lines change the graph colour over time
    #     # each RGB channel follows a different sinusoid, with the peaks of all
    #     #   three seperated to provide periodic chromatic shifting
    #     # the amplitude is in order to prevent unaesthethic light colours
    #     r = abs(0.75*np.cos(0.5 * (i / 100 * 2 * np.pi + np.pi)))     
    #     g = abs(0.75*np.sin(0.5 * (i / 100 * 2 * np.pi) + 0.2 * np.pi))
    #     b = 0.7
    #     # plot oscillators' displacement at each time step
    #     # use circles on string to represent the N oscillators
    #     ax.plot(positions, slices[i], marker='o', markersize=3, color=(r, g, b), label=f'time t = {str(i) + "ms" if i < 1000 else str(round(i / 1000, 3)) + "s"}')
    #     ax.legend(loc="upper right")

    #     ax.set_xlim([0,L])
    #     ax.set_ylim([-amp_0, amp_0])

    # ani = FuncAnimation(fig, animate, frames = n_step, interval = 1, repeat=False)
    # annotate graph and display it to user
    
    # plt.title(f"Oscillation of lattice over time with non-linear coefficient alpha = {alpha}")
    # plt.xlabel("Position along string, p [m]")
    # plt.ylabel("Horizontal displacement, x [m]")

    # plt.show()


f = open("userdata.txt", "r")
data = f.readlines()
for i in range(len(data)):
    data[i] = data[i].strip()

Fermi({'L': data[0], 'tf': data[1], 'N': data[2], 'k': data[3], 'rho': data[4], 'alpha': data[5], 'amp_0': data[6], 'shape_0': data[7]})