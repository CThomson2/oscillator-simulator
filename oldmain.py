import math
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
from scipy.integrate import solve_ivp
import webbrowser
import os

# import modules
from methods import Methods as mt

# Boundary Conditions
# since x_0(t) = x_L(t) = 0, velocity at the boundaries is also a constant zero
x0 = 0
xL = 0
v0 = 0
vL = 0

class String():

    # ----- Parameter Definitions -----
    # x_i: displacement from equilibrium of i-th oscillator
    # u_i: velocity of i-th oscillator (eqbm vel. is zero)
    # p_i: distance along string (can be represented as i*h)
    # t: time
    # h: distance between oscillators (step size)
    # L: length of string
    # N: number of oscillators

    def __init__(self, data):
        self.data = data
        
        self.L, self.tf, self.N, self.k, self.rho, self.alpha, self.pert_ini = self.destruct_dict(self.data, 'L', 'tf', 'N', 'k', 'rho', 'alpha', 'pert_ini')
        self.N = int(self.N)
        self.tf = int(self.tf)

        # calculate stepsize (distance between oscillators) and mass of each oscillator
        self.h = self.L / (self.N - 1)  # step size
        self.m = self.rho * self.h  # mass of osc.
        # array holding positions of oscillators along string
        # will be used as the x-axis for plotting variables of interest along length of string
        self.pos_lattice = [i * self.h for i in range(self.N)]
        # choose a timestep (also in ms)

        # IMPORTANT - dt shouldn't be hardset to 1. Sometimes that will be too large.
        self.dt = 0.01
        # But, if it's less than 1, the t array will contain decimals (e.g. t = [0, 0.2, 0.4, ... , 99.8, 100])
        # Clearly, you'll need to normalise the times so that the time step is equal to 1, which will increase
        # the final time. This means that when you plot the time, you need to use the original timescale provided
        # by the user.

        self.n_step = math.ceil(self.tf / self.dt)
        # create time-step array
        self.t = np.linspace(0, self.tf - 1, self.n_step, dtype=int)
        # create arrays that will hold all displacements, velocities and accelerations in the 2D space of discretised time and distance along string
        self.x = []
        self.v = []
        # 2D array will hold an array of displacements for each time-step
        for i in range(self.n_step):
            self.x.append(np.zeros(self.N))
            self.v.append(np.zeros(self.N))

    # Create a class method that breaks the dictionary of user input down into its constituent parameters
    # Will be used to extract relevant parameters for each method
    @classmethod
    def destruct_dict(cls, dict, *args):
        return (float(dict[arg]) for arg in args)

    # -----
    # Aim: to find the velocity of each oscillator at each point in time
    # Method: create a 2D array to store arrays of velocities at each time-step for each oscillator

    # STEP 1: create function to find initial displacement from user-specified conditions
    # Additionally, store initial velocities and accelerations in main solution array
    # -----

    def get_init(self):
        init_type = self.data['init_type']

        x_init = np.zeros(N)
        # --- Case 1 ---
        if init_type == 'sine':
            for i in range(N):

                x_init[i] = self.pert_ini * math.sin((2 * math.pi / L) * self.pos_lattice[i])
                # handle the potential error of boundary values being <just off> zero when they should be zero exactly
                if abs(x_init[i]) < 1e-10:
                    x_init[i] = 0
            return x_init
        # --- Case 2 ---
        elif init_type == 'half-sine':
            for i in range(N):

                x_init[i] = self.pert_ini * math.sin((math.pi / L) * self.pos_lattice[i])
                if abs(x_init[i]) < 1e-10:
                    x_init[i] = 0
            return x_init
        # --- Case 3 ---
        elif init_type == 'parabola':
            for i in range(N):

                x_init[i] = -(4 * self.pert_ini / L ** 2) * \
                    (i * self.h - L / 2) ** 2 + self.pert_ini
                if abs(x_init[i]) < 1e-10:
                    x_init[i] = 0
            return x_init
        # Add more polynomial possibilites here (cubic, quartic etc.)
    

    # Create method to store initial conditions in first time-step of arrays
    def setup(self):
        init_disp = self.get_init()

        plt.plot([i * self.h for i in range(N)], init_disp, 'r.-')
        plt.savefig('./static/img/graph.png')
        plt.cla()


        # Set the oscillators' displacement at time t=0 to initial conditions
        self.x[0] = init_disp

        # # Clearly, the initial velocity is zero at every point
        init_vel = np.zeros(self.N)
        self.v[0] = init_vel


    # ----
    # STEP 2: Define governing equations and derivative functions
    # ----

    # define governing acceleration equation
    def get_acc(self, osc):
        result = (self.k / self.h) / self.m * (osc[2] + osc[0] - 2 * osc[1]) * (1 + self.alpha * (osc[2] - osc[0]))
        return result

    def model(self, t, y):
        x = y[:int(len(y) / 2)]
        u = y[int(len(y) / 2):]

        fx = u
        fu = np.zeros(N)

        for i in range(1, N-1):
            fu[i] = self.get_acc([ x[i - 1], x[i], x[i + 1] ])

        return np.concatenate([fx, fu])

    
    # -----
    # STEP 3: Taking slices of time, we study the velocities and displacements of each oscillator using provided recurrence relations
    # -----
    def iterate(self):
    # -----

        """To find displacement at any given point, we need its derivative: velocity.
        To find velocity, we need the acceleration at that point.
        Thus, we need to fill the velocity and displacement arrays iteratively.
        i.e. we find velocity at point i at time-step j using the formula for acceleration at point i at time-step j-1
        We do this for all oscillators i in-between the boundary conditions, which remain at zero displacement throughout
        Next, we find displacement at point i at time-step j using the indexed value for velocity that we just found and
        stored in the appropriate location in the global velocity array.
        Before this algorithm, we can fill some rows in the 2D grids (arrays) for v and x
        We have already inserted the initial conditions, but we also know the displacement at time 0 + ∂t (i.e. first time step)
        Since v_0(t) = 0, the displacement cannot change in the second time-step. Therefore x_i(t1) = x_i(t0)"""

        
        # Initial conditions
        t0 = 0
        y0 = np.concatenate([self.x[0], self.v[0]])

        # solve simt. equations using solve_ivp
        t_eval = np.linspace(0, self.tf, num=3000)
        
        res = solve_ivp(self.model, [0, self.tf], y0)

        slices = []
        for i in range(len(res.y[0])):
            slices.append(np.zeros(self.N))
            for j in range(self.N):
                slices[i][j] = res.y[j][i]

        return slices


    def plot(self, data):
        fig, ax = plt.subplots()

        print(len(data))

        def animate(i):
            print(i)
            ax.clear()

            r = abs(np.cos(0.5 * (i / 100 * 2 * np.pi + np.pi)))
            g = abs(np.sin(0.5 * (i / 100 * 2 * np.pi) + 0.2 * np.pi))
            b = abs(np.sin(0.5 * (i / 100 * 2 * np.pi - 0.75 * np.pi)))
            ax.plot(self.pos_lattice, data[i], color=(r, g, b))

            ax.set_xlim([0,L])
            ax.set_ylim([-self.pert_ini, self.pert_ini])
        
        ani = FuncAnimation(fig, animate, frames = len(data), interval = 25, repeat=False)

        plt.xlabel("Position along lattice")
        plt.ylabel("Oscillator displacement")
        plt.show()

    def run(self):
        self.get_init()
        self.setup()
        data = self.iterate()
        self.plot(data)


# Hardcoded parameters for testing and debugging
L = 1
t_final = 100
N = 100
k = 0.1
rho = 400
alpha = 0.5
pert_ini = 1.0
init_type = 'sine'

data = {'L': L, 'tf': t_final, 'N': N, 'k': k, 'rho': rho,
         'alpha': alpha, 'pert_ini': pert_ini, 'init_type': init_type}

mystr = String(data)
mystr.run()