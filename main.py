# import packages
import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
import webbrowser
import os

# import modules
from methods import Methods as mt

print(mt.euler_simple())

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
        self.dt = 1
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
        self.a = []
        # 2D array will hold an array of displacements for each time-step
        for i in range(self.n_step):
            self.x.append(np.zeros(self.N))
            self.v.append(np.zeros(self.N))
            self.a.append(np.zeros(self.N))

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

        # -----
        # show user the initial conditions in graphical form
        # filename = 'file:///' + os.getcwd() + '/' + 'initcond.html'
        # webbrowser.open_new_tab(filename)
        # -----

        # Set the oscillators' displacement at time t=0 to initial conditions
        self.x[0] = init_disp

        # # Clearly, the initial velocity is zero at every point
        init_vel = np.zeros(self.N)
        self.v[0] = init_vel

        # # Initial acceleration can be found using values for displacement
        self.a[0][0] = 0
        self.a[0][self.N - 1] = 0
        for i in range(1, self.N - 1):
            self.a[0][i] = self.get_acc([self.x[0][i - 1], self.x[0][i], self.x[0][i + 1]])

        # plt.figure(figsize=(6, 3))
        # plt.subplot(131)
        # plt.plot(self.pos_lattice, self.x[0], 'b.-')
        # plt.xlabel('string')
        # plt.ylabel('displacement')
        # plt.subplot(132)
        # plt.plot(self.pos_lattice, self.a[0], 'r.-')
        # plt.xlabel('string')
        # plt.ylabel('acceleration')
        # plt.show()

    # print('Mass:', m, 'Step size h:', h, 'Length density rho:', rho)
    # print(
    #     f'Init. Disp: {x[0]}, \n\nInit. Vel: {v[0]}, \n\nInit. Acc: {a[0]}\n')

    # ----
    # STEP 2: Define governing equations and derivative functions
    # ----

    # define governing acceleration equation
    def get_acc(self, osc):

        # check for boundary conditions
        # for ox in osc:
        #     if ox in [0, self.N - 1]:
        #         return 0
        # otherwise use i parameter to determine the acceleration of the i-th oscillator

        result = (self.k / self.h) / self.m * (osc[2] + osc[0] - 2 * osc[1]) * (1 + self.alpha * (osc[2] - osc[0]))
        return result

    def model(self, t, y):
        # print(y[1:3], y[-4:-2])
        x = y[:int(len(y) / 2)]
        v = y[int(len(y) / 2):]
        fx = v
        fv = np.zeros(len(v))

        for i in range(len(v)):
            if i not in [0, len(v) - 1]:
                fv[i] = self.get_acc(
                    [x[i - 1], x[i], x[i + 1]]
                )
            elif i == 0:
                fv[i] = self.get_acc(
                    [x0, x[i], x[i + 1]]
                )
            elif i == len(v) - 1:
                fv[i] = self.get_acc(
                    [x[i - 1], x[i], xL]
                )
        # print(fv)
        return np.concatenate([fx, fv])

    
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

    # print(f'\n{t}\n')

    # throughout this program, i is associated with position along string and j with time-step
    # pattern is continued here for clarity and consistency
        
        # Initial conditions
        t0 = 0
        y0 = np.concatenate([self.x[0][1:-1], self.v[0][1:-1]])

        print(len(y0))
        # print(y0)

        # solve simt. equations using solve_ivp
        t_eval = np.linspace(0, self.tf, num=3000)
        res = solve_ivp(self.model, [t0, self.tf], y0, t_eval=t_eval)

        # print('\n length of res.y:', len(res.y))

        # create list of slices in time wherein we store all displacements at that instant
        timesteps = []
        for i in range(0, len(res.y[0])):
            # store each instant in time in timesteps array
            timesteps.append([j[i] for j in res.y[int(len(res.y) / 2):]]) # dividing by 2 to discard velocities

        # plot displacement for first 50 timesteps
        for i in range(50):
            plt.plot( self.pos_lattice, np.concatenate([[x0], timesteps[i], [xL]]) )
        plt.title(f't_eval parameter = {len(t_eval)}, alpha = {self.alpha}')
        plt.show()
        
        # for j in range(0, self.tf, 10):
        #     plt.plot(self.pos_lattice, np.concatenate([[x0], timesteps[j], [xL]]))
        #     plt.show()

    def iterate_old():
        # iterating through time-steps
        # for j in range(0, self.n_step - 1):
            # iterating through oscillators on string
    
            # Euler's Method
            # for i in range(1, N - 1):
            #     v_slope = self.get_acc([self.x[j][i - 1], self.x[j][i], self.x[j][i + 1]])
            #     self.v[j + 1][i] = self.v[j][i] + v_slope * self.dt
            #     self.x[j + 1][i] = self.x[j][i] + self.v[j][i] * self.dt

            # SimEqns Method
            # f = model
            

            # # Runge Kutta Method
            # for i in range(1, N - 1):
            #     # compute slope using the ODE
            #     t_d = self.t[j]
            #     v_d = self.v[j][i]
            #     k1 = self.get_acc(j, i)

            #     t_d = self.t[j] + self.dt / 2
            #     v_d = self.v[j][i] + k1 * self.h / 2
            #     k2 = self.get_acc(j, i)

            #     t_d = self.t[j]
            #     v_d = self.v[j][i]
            #     k1 = self.get_acc(j, i)

            #     t_d = self.t[j]
            #     v_d = self.v[j][i]
            #     k1 = self.get_acc(j, i)
                

            #     x_dummy = x_rk[i] + h / 2
            #     y_dummy = y_rk[i] + k1 * h / 2
            #     k2 = model(y_dummy, x_dummy)

            #     x_dummy = x_rk[i] + h / 2
            #     y_dummy = y_rk[i] + k2 * h / 2
            #     k3 = model(y_dummy, x_dummy)

            #     x_dummy = x_rk[i] + h
            #     y_dummy = y_rk[i] + k3 * h
            #     k4 = model(y_dummy, x_dummy)

            #     # compute slope as weighted average of four slopes
            #     slope = 1 / 6 * k1 + 2 / 6 * k2 + 2 / 6 * k3 + 1 / 6 * k4

            #     # use the RK method
            #     y_rk[i + 1] = y_rk[i] + h * slope

            # Fv = lambda v_next, i: self.v[j][i] + self.get_acc(self.t[j + 1], i) * self.dt - v_next
            # # print('\n'+str(j)+'\n')
            # for i in range(1, N - 1):
            #     self.v[j + 1][i] = mt.muller(Fv, self.v[j][i], 1.01 * self.v[j][i] + 10 ** -3, 0.98 * self.v[j][i], i, 100)
            #     self.x[j + 1][i] = self.x[j][i] + self.v[j][i] * self.dt
        pass


    # Plot results
    def plot(self):
        # b1 = {
        #     "x": [],
        #     "v": []
        # }
        # for ball in x:
        #     b1["x"].append(ball[4])
        # for ball in v:
        #     b1["v"].append(ball[4])

        # plt.plot(t, b1["x"], 'r.-', t, b1["v"], 'g.-')

        # plt.figure(figsize=(9, 6))

        # plt.subplot(131)
        # plt.plot(pos_lattice, x[0])
        # plt.subplot(132)
        # plt.axes((0, -1.0, 1.0, 2.0))
        # plt.plot(pos_lattice, x[100])
        # # plt.axes((0, -1.0, 1.0, 1.0))
        # plt.subplot(133)
        # plt.plot(pos_lattice, x[200])
        # plt.axes((0, -1.0, 1.0, 1.0))
        # plt.subplot(1, 3, (2, 1))
        # plt.plot(pos_lattice, x[3])
        # plt.subplot(1, 3, (2, 2))
        # plt.plot(pos_lattice, x[4])
        pass

    def run(self):
        self.get_init()
        self.setup()
        self.iterate()
        
        # self.plot()


# Hardcoded parameters for testing and debugging
L = 1
t_final = 10000
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