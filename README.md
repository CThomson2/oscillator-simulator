# FermiPastaG12

Simulation project for an engineering numerical methods class. The program will solve a dual system of ODEs and find displacements of oscillators in a 1D lattice with user-specified parameters, display an animation of the oscillators' displacements through time and finally perform a Fourier analysis.

__## [Video Demonstration (Non-Linear Case)](https://vimeo.com/768030094)__

# Simulation Program Structure

## HTML GUI

### Two HTML forms

#### Form 1

- Has every input apart from initial conditions, x(t0) for all i.
- Hides submit button until every input field contains a number.

#### Form 2

- First asks what shape the string will take initially. Options are limited to sine, half-sine or parabola (since x0 and xL must be zero).
- If user inputs "sine" or "half-sine", ask for amplitude and number of wavelengths or half-wavelengths, resp. Check no. of wavelengths is an integer.
- If user inputs "parabola", ask for stationary point height.
  - Extra feature: _Add another option: "polynomial"_:
    - Allow user to chooose polynomial degree and coefficients
      - i.e. if degree is "4", ask user for five coefficients (in order)
    - allow user to click button to randomize coefficients before/after choosing polynomial degree
    - keep button shaded and unclickable until polynomial satisfies the constraints of x0 = xL = 0
      - or, if polynomial doesn't fit constraints, translate it until it does without changing (challenging)

#### Python integration

- Take data from user and call main "simulate()" function with paramaters
- Represent processed data and solution on new HTML page to show user

---

## Solve System of ODEs

### Simultaneous Equations

- dx_i/dt for every oscillator i, i.e. N simultaneous equations

### Perform RungeKutta OR Euler method on d^2x/dt^2 to determine discretised velocities

- #### Create array of arrays, where each array holds all velocity values for one of N oscillators

### Repeat process for dx/dt for displacements

- #### Instead of using a derivative function, we'll have to index the velocity arrrays to find dx_i/dt at each time-step

---

## Plot x(t, i)

- X-axis will be the distance along string (i-value)
- To show the progression through time, we will either animate the graph or plot M graphs in one figure, where M is the number of time-steps between t_0 and t_final
  - Matplotlib has an animation class

---

## Analyse the change in the first few Fourier Series coefficients over time in the non-linear case

---

## Additional Features

- If user inputs a step size large enough to cause instability, kill the program and print error OR even prevent user from submitting the Parameters form if their h value is too high - we will need to work out k first to determine the region of stability
- Create a range of preset polynomials that satisfy the boundary conditions for the user to choose from when setting initial conditions. Create a button that when clicked unhides a subsection of polynomial graphs that can be clicked.
