# Import modules
from re import L
import tqdm
import argparse
import math
import random
import numpy as np
import matplotlib.pyplot as plt
import functions as func

# Parameter definition (default)
Temp_initial = 1000
Temp_min = 500
Markov_length = 60 #Internal calculation
damping_factor = 0.6

num_units = 10
num_feat  = 10
max_value = 1
min_value = 0

# Function definition
def objective_fun(sel_flag):
    """This function will return the chosen objective"""
    # Select the function 
    if (sel_flag == "ackley"):
        return func.ackley
    elif (sel_flag == "beale"):
        return func.beale
    elif (sel_flag == "booth"):
        return func.booth
    elif (sel_flag == "bukin6"):
        return func.bukin6
    elif (sel_flag == "crossintray"):
        return func.crossintray
    elif (sel_flag == "easom"):
        return func.easom
    elif (sel_flag == "eggholder"):
        return func.eggholder
    elif (sel_flag == "goldstein"):
        return func.goldstein
    elif (sel_flag == "himmelblau"):
        return func.himmelblau
    elif (sel_flag == "holdertable"):
        return func.holdertable
    elif (sel_flag == "levi"):
        return func.levi
    elif (sel_flag == "matyas"):
        return func.matyas
    elif (sel_flag == "rastringin"):
        return func.rastrigin
    elif (sel_flag == "rosenbrock"):
        return func.rosenbrock
    elif (sel_flag == "schaffer2"):
        return func.schaffer2
    elif (sel_flag == "sphere"):
        return func.sphere
    elif (sel_flag == "threehump"):
        return func.threehump
    
def matropolis(diff, temp):
    """Check whether we should accept the result or not """
    if diff < 0:
        return 1
    else:
        prob = math.exp(-diff / temp)
        if prob > random.rand():
            return 1
        else:
            return 0

def gen_new(x_old, dis):
    """ this function will generate new value based on the current
    temperature and the type we want to add.
    ---
    Args:
        x_old (n_particles, dimensions): old x
        dis (int): used for updating the position
    """
    
    x_new = np.copy(x_old)

    for i in range(x_new.shape[1]):
        x_new[0][i] = x_new[0][i] + np.random.uniform(low=-0.055, high=0.55) * dis

    return x_new
            
    
def sa_simulation(sel_flag, distance,x_init, temp_init = Temp_initial, temp_min = Temp_min, markov_len = Markov_length, damping_factor = damping_factor, max_value = max_value, min_value = min_value):
    """ sa_simulation function: 
    This function will do the simulated annealing simulation
    -------
    Args:
        sel_flag (string): choose the objective function
        distance (int) : used for generating new positions
        x_init (num_dim, num_feat) : initial value of x
        temp_init (int): Initial temperature
        temp_min (int): End temperature
        markov_len (int): Internal calculation counts
        damping_factor (float): rate of decreasing temperature
    """
    # Param definition
    Temp = temp_init
    x_ = np.atleast_2d(x_init)
    
    # Random generator
    randseed = random.randint(1,50)
    random.seed(randseed)
    
    # Sel objective function
    obj_fun = objective_fun(sel_flag)
    
    # Run the simulation
    while Temp >temp_min:
        for i in range(markov_len):
            # Calculate fitness
            y_old = obj_fun(x_)
            # Generate new x
            x_new = gen_new(x_, distance)
            y_new = obj_fun(x_new)
            
            # Selection
            if (matropolis(y_new - y_old, Temp)):
                x_ = x_new
        
        # Update the temp
        Temp = Temp * damping_factor
        
    return x_
    
    