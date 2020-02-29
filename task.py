""" CREATE YOUR SOLVER

Use tutorial.py as a reference

SOLVER
1. Import the correct modules
    a) solver, FileLoader and DataCollector
2. Input the initial variables to create an instance of the solver
    a) make the length of the solver 3 metres (inputting 3 will work)
    b) calculate and set the wavelength using omega_0 (remember that c = frequency * wavelength)
    c) make the simulation time = size / C (think about what this means, the speed of the wave is C)
    d) make eps_r_max = 3
3. Create a pulse using add_oscillating_pulse, store it in a variable called pulse (replace None)
4. Plot the pulse (run the code)
5. Create your material and store it in a variable called material using create_material (replace None)
6. Set the upper boundary to be absorbing using set_reflect_boundaries
7. Add a square to your material in the bottom right, upper left corner (2.4, 2.4), lower right corner (2.9, 2.9)
8. Save the file as file name 'task'
9. Call the solve function (finally), setting realtime to be False
10. When the loading has reached 100%, check the directory that this file is located in, you should see a file
   called task.npy . After this file has been created, you no longer need to run or use solver and can move on to
   using the File Loader and Data Collector for analysis

FILE LOADER
1. Comment out the code relating to the solver
2. Create a variable called file_loader and create an instance of the FileLoader, giving the correct file name
3. Use the play function of file loader, setting interval=50
    a) Change the value of interval from 1 to 200, see the effects of the animation

DATA COLLECTOR
1. coming soon...

"""

# PART 1

# initiate variables
sigma_w = 0.5 * 10 ** 9  # frequency bandwidth, this is part of the GAUSSIAN envelope
omega_0 = 10 * 10 ** 9  # central frequency, this is how fast your pulse will oscillate
s = 10  # mesh points per wavelength
stability = 0.1  # time mesh stability factor
C = 2.99 * (10 ** 8)


# PART 2
wavelength = None
size = None
time = None


# PART 3
pulse = None
# PART 4

# PART 5
material = None

# PART 6, 7, 8, 9, 10