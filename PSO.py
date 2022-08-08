import argparse
import tqdm
import numpy as np
import matplotlib.pyplot as plt
from functions import *

# the function that needs to be minimised


def function(x):
    return np.sum(np.square(x - 2), axis=1, keepdims=True)


# hyperparameters
def main():
    # initialisation
    args = _parse_args()
    num_units = args.num_units
    num_feat = args.num_feat
    position = np.random.random((num_units, num_feat))
    """

    #Allows the numbers to start from a specific range
    for i in range(num_feat):
        a=float(input("Enter the maximum value "+str(i+1) +" point can have:"))
        b=float(input("Enter the minimum value "+str(i+1) +" point can have:"))
        position[:,i]=position[:,i]*(a-b)+b
    """
    best_individual = np.copy(position)
    best_pos = position[np.argmin(function(best_individual)), :].reshape(1, num_feat)
    velocity = np.random.random((num_units, num_feat))
    valuesnow = function(position)
    indibestnow = np.copy(valuesnow)
    # Iteration to reach the stable point
    for i in tqdm.tqdm(range(args.iters)):
        r1 = np.random.random()
        r2 = np.random.random()
        velocity = (
            r1 * velocity
            - args.c1 * r1 * (position - best_individual)
            - args.c2 * r2 * (position - best_pos)
        )
        position = position + velocity
        valuesnow = function(position)
        temp = valuesnow < indibestnow
        temp = temp.reshape(num_units, 1)  # for the eggholder
        best_individual = best_individual + temp * (position - best_individual)
        indibestnow = indibestnow + temp * (valuesnow - indibestnow)
        best_pos = best_individual[np.argmin(function(best_individual)), :].reshape(
            1, num_feat
        )
        # If it reaches a stable point which is not the maxima, it will cause an unequilibruim
        if np.isclose(np.max(valuesnow), np.min(valuesnow), rtol=1e-4, atol=1e-8):
            velocity = np.random.random((num_units, num_feat)) * best_pos
        # Getting a visual for how particles move for two dimensions
        # valuesnow=function(position)
    print(best_pos)
    print(np.min(indibestnow))


def _parse_args():
    """
    The different arguments for the parser
    """
    parser = argparse.ArgumentParser(description="Obtain the parameters for the PSO")
    parser.add_argument(
        "-d",
        "--num_feat",
        type=int,
        default=10,
        help="Number of dimensions each particle has",
    )
    parser.add_argument(
        "-n",
        "--num_units",
        type=int,
        default=10,
        help="Number of particles for the run",
    )
    parser.add_argument("--c1", type=float, default=2, help="Positional Dependence")
    parser.add_argument("--c2", type=float, default=2, help="Global Dependence")
    parser.add_argument("--modified", action="store_true", help="Use the modified PSO")
    parser.add_argument(
        "--iters",
        type=int,
        default=1000,
        help="Number of iterations the algorithm will run",
    )
    parser.add_argument(
        "--func", type=str, default="sqrt_loss", choices=["sqrt_loss", "eggholder_loss"]
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
