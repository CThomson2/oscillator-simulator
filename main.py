# import packages
import math
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
from scipy.integrate import solve_ivp
from scipy.fft import fft
from scipy.signal import find_peaks
from scipy import interpolate
from random import randint

# import modules
from initial_conditions import get_init

###

# ----- Introduction -----

### Program Capabilities and Limitations

# will simulate 1D FPUT lattice for any numerical combination of initial conditions,
#   however in its current version only allows for initial waveforms of the sine,
#   half-sine and parabolic waves
# final animation takes approx. 10-15 minutes to render fully (tested on an 8GB 2021 iMac),
#   therefore the full Fourier coefficient graph will be displayed in the written report
# contains an algorithm that allows position and velocity calculation of only half of the N
#   oscillators (sine wave only), as the other half are simply the negative of these values

###

# define colour scheme for modal energy graph, to be used when plotting
MPL_COLOUR_SCHEME = {
    "Mode1": (42/255, 157/255, 143/255),
    "Mode2": (233/255, 196/255, 106/255),
    "Mode3": (244/255, 162/255, 97/255),
    "Mode4": (231/255, 111/255, 81/255)
}

def fermi_sim(data):

    # ----- DEFINE PARAMETERS, USER DATA AND INITIAL CONDITIONS ----- #
    # --------------------------------------------------------------- #

    # create a function that breaks the dictionary of user inputs down into its constituent parameters
    destruct_dict = lambda dict, *args: (float(dict[arg]) for arg in args)

    # extract all variables from the data dictionary
    L, tf, N, k, rho, alpha, amp_0 = destruct_dict(data, 'L', 'tf', 'N', 'k', 'rho', 'alpha', 'amp_0')
    shape_0 = data['shape_0']
    anim = data['anim']

    # set following params to integer format so they can be iterated through later
    N = int(N)
    tf = int(tf)

    # calculate total time steps the variables of interest will be recorded at
    dt = 1
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


    # ----- RUNGE KUTTA LOOP ----- #
    # ---------------------------- #

    # define model function to be used as an argument in solve_ivp later
    # this function returns the derivatives (slopes) of all input variables
    def model(t, y):
        
        # to minimise the computations required, we make use of the sine wave's symmetry - present
        #   even in the non-linear case!
        if shape_0 == 'sine':

            # instead of evaluating the slopes of all N oscillators, we only need to do so for the first
            #   half, then invert the values for the second half of the wave
            # there is a notable difference in the algorithm depending on N's parity (i.e. odd or even)
            if N % 2 == 0:
                # extract the displacements and velocities from y array
                x = y[:int(len(y) / 4)]
                u = y[int(len(y) / 2):int(len(y) * 3 / 4)]

                # the derivative of displacement is simply the velocity
                fx = u
                fu = np.zeros(math.floor(N/2))

                # fill dv/dt array using acceleration function and provided displacement values
                # begin range at index 1 because of the boundary condition
                for i in range(1, len(fu) - 1):
                    fu[i] = get_acc([ x[i - 1], x[i], x[i + 1] ])
                # -x[-1] is simply the next oscillator in line, and lies on the other side of the x-axis
                #   as x[-1]
                fu[-1] = get_acc([ x[-2], x[-1], -x[-1] ])

                # concatenate the derivaties to the the corresponding reversed negative derivatives (envision
                #   the sine wave if this is unclear)
                return np.concatenate([fx, [-f for f in reversed(fx)], fu, [-f for f in reversed(fu)]])
            
            # when N is odd, the central oscillator remains static throughout
            # therefore, to find the acceleration value for the final velocity value
            #   we must simply use zero as the displacement of the subsequent oscillator
            else:
                x = y[:int(len(y) / 4)]
                u = y[int(len(y) / 2):int(len(y) * 3 / 4)]

                fx = u
                fu = np.zeros(math.floor(N/2))

                for i in range(1, len(fu) - 1):
                    fu[i] = get_acc([ x[i - 1], x[i], x[i + 1] ])
                fu[-1] = get_acc([ x[-2], x[-1], 0 ])

                return np.concatenate([fx, [0], [-f for f in reversed(fx)], fu, [0], [-f for f in reversed(fu)]])
        
        # unfortunately, this mathematical exploit is invalid for the non-linear half-sine or parabolic cases
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


    # ----- FOURIER ANALYSIS ----- #
    # ---------------------------- #

    # create array of fourier coefficients for the N oscillators at every time step
    fourier = []
    for i in range(0, len(slices)):
        # since the coefficients represent the amplitude of each frequency,
        #   the magnitude is all that we care about.
        # normalise coefficients by dividing by half of N.
        fourier.append(np.absolute(np.fft.rfft(slices[i])) / (N / 2))

    # we now have the Fourier series for every timestep, but we want to display the evolution
    #   of the coefficients over time.
    # invert the Fourier series by swapping i-th, j-th values with their converse and appending
    #   to a new coefficients array
    coefs = []
    for j in range(len(fourier[1])):
        coefs.append(np.zeros(len(fourier)))
        for i in range(len(fourier)):
            coefs[j][i] = fourier[i][j]
  
    # find the timestamps at which point the waveform peaks
    # create a list of these timestamps for all four coefficients
    peak_energies = []
    # iterating from 1 to 5 as we aren't plotting as the first Fourier Coefficient is invalid (zero)
    for i in range(1, 5):
        # use scipy.signal's method for finding local peaks; this returns the time values of the peaks
        peaks = list(find_peaks(coefs[i])[0])
        peak_energies.append(peaks)


    # ----- GRAPHING SOLUTION ----- #
    # ----------------------------- #

    # the final stage is plotting the results
    # if the user selected the animation version, it is achieved by using MatPlotLib's animation class,
    #   where the solution is plotted on a live graph
    #   otherwise, plot the static Fourier graph and a figure showing the solution at seperate snapshots in time

    if not anim:
        # plot Fourier coefficients (modal energy over time)
        plt.plot(
            peak_energies[0], [coefs[1][i] for i in peak_energies[0]],
            color=MPL_COLOUR_SCHEME["Mode1"], label='Mode 1'
        )
        plt.plot(
            peak_energies[1], [coefs[2][i] for i in peak_energies[1]],
            color=MPL_COLOUR_SCHEME["Mode2"], label='Mode 2'
        )
        plt.plot(
            peak_energies[2], [coefs[3][i] for i in peak_energies[2]],
            color=MPL_COLOUR_SCHEME["Mode3"], label='Mode 3'
        )
        plt.plot(
            peak_energies[3], [coefs[4][i] for i in peak_energies[3]],
            color=MPL_COLOUR_SCHEME["Mode4"], label='Mode 4'
        )
        plt.legend(loc="upper center")

        plt.xlim([0, tf])
        plt.ylim([0, 1])

        plt.title("Fourier series of lattice oscillations at discrete time steps")
        plt.xlabel("time, t [ms]")
        plt.ylabel("Fourier Coefficients")

        # plot series of time-seperated lattice displacements
        fig = plt.figure(figsize=(15,12))
        # create list of timestamps
        snapshots = [slices[int(t)] for t in np.linspace(0, 300000-15000, 20)]
        tstep = 0
        for i in range(1, 21):
            plt.subplot(4, 5, (i, i))
            plt.plot(positions, snapshots[tstep], marker='o', markersize=3, color=MPL_COLOUR_SCHEME["Mode4"],
                label=f'time t = {str(tstep * 600) + "ms" if tstep <= 1 else str(round(tstep * 600 / 1000, 3)) + "s"}')
            tstep += 1
            plt.ylim([-amp_0, amp_0])
            plt.legend(loc="upper left")
        
        plt.show()
        
        fig.suptitle("Full time evolution of lattice with non-linear term alpha =", alpha)
        fig.tight_layout(pad=3.0)
        # once we've plotted the above figures, the simulation is complete
        return

    # set up MPL figure with two axes (f: fourier and s: solution)
    fig, (ax_f, ax_s) = plt.subplots(2, 1)
    fig.set_size_inches(12, 10)

    # using a doubele list comprehension, create a multi-dimensional list that stores,
    #   for each coefficient, the peak f(t) values
    c = [[coefs[i + 1][j] for j in peak_energies[i]] for i in range(0, 4)]

    # create spline arrays for all coefficients such that they equal n_step in length and can 
    #   thus be animated over for n_step frames
    c_splines = []
    for i in range(len(c)):
        f_spline = interpolate.splrep(peak_energies[i], c[i], s=0)
        c_splines.append(interpolate.splev(t_eval, f_spline, der=0))


    def animate(i):
        # reset the graphs at each step
        ax_s.clear()
        
        # --- x_i(t) Solution ---

        # the following lines change the graph colour over time
        # each RGB channel follows a different sinusoid, with the peaks of all
        #   three seperated to provide periodic chromatic shifting
        # the 0.75 amplitude is in order to prevent unaesthethic light colours
        r = abs(0.75*np.cos(0.5 * (i / 100 * 2 * np.pi + np.pi)))     
        g = abs(0.75*np.sin(0.5 * (i / 100 * 2 * np.pi) + 0.2 * np.pi))
        b = 0.7
        # plot oscillators' displacement at each time step
        # use circles on string to represent the N oscillators
        ax_s.plot(positions, slices[i], marker='o', markersize=3, color=(r, g, b), label=f'time t = {str(i) + "ms" if i < 1000 else str(round(i / 1000, 3)) + "s"}')
        ax_s.legend(loc="upper right")
        ax_s.set_xlim([0,L])
        ax_s.set_ylim([-amp_0, amp_0])
        plt.fill_between(positions, -1, slices[i], color=MPL_COLOUR_SCHEME["Mode4"])


        # annotate figure
        ax_s.set_title(f"Oscillation of lattice over time with non-linear coefficient alpha = {alpha}")
        ax_s.set_xlabel("Position along string, p [m]")
        ax_s.set_ylabel("Horizontal displacement, x [m]")

        # --- Fourier Coefficients ---
        ax_f.clear()

        # plot the coefficients of interest
        ax_f.plot(c_splines[0][:i], color=MPL_COLOUR_SCHEME["Mode1"], label='Mode 1')
        ax_f.plot(c_splines[1][:i], color=MPL_COLOUR_SCHEME["Mode2"], label='Mode 2')
        ax_f.plot(c_splines[2][:i], color=MPL_COLOUR_SCHEME["Mode3"], label='Mode 3')
        ax_f.plot(c_splines[3][:i], color=MPL_COLOUR_SCHEME["Mode4"], label='Mode 4')
        ax_f.legend()

        ax_f.set_xlim([0, tf])
        ax_f.set_ylim([0, 1])

        ax_f.set_title("Fourier series of lattice oscillations at discrete time steps")
        ax_f.set_xlabel("time, t [ms]")
        ax_f.set_ylabel("Fourier Coefficients")
            

    ani = FuncAnimation(fig, animate, frames = n_step, interval = 1, repeat=False)

    # display figure to user
    fig.tight_layout(pad=5.0)
    plt.show()



# firstly, the program reads the data from the file view.py, where we submit the user's inputs to
# this data is saved to a list that is used as the parameter to the main simulation function
f = open("userdata.txt", "r")
data = f.readlines()
for i in range(len(data)):
    data[i] = data[i].strip()

# change the "checked" value of the animate option checkbox to a Boolean for easier manipulation
# if checkbox was checked, it will be recorded as "on". If not, it won't be recorded and data will
#   only have eight items
if len(data) == 9:
    data[8] = True
else:
    data.append(False)

fermi_sim({'L': data[0], 'tf': data[1], 'N': data[2], 'k': data[3], 'rho': data[4], 'alpha': data[5], 'amp_0': data[6], 'shape_0': data[7], 'anim': data[8]})