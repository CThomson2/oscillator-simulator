# Simulation Program Structure

## HTML GUI

### Two HTML forms.

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

## Solve System of ODEs
