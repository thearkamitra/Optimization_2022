import argparse
import tqdm
import numpy as np
import matplotlib.pyplot as plt
import functions as func

# the function that needs to be minimised


def f(x):
    return np.sum(np.square(x - 2), axis=1, keepdims=True)


# hyperparameters
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
    parser.add_argument(
        "--max_value",
        type=float,
        default=1,
        help="Maximum value of the elements initially",
    )
    parser.add_argument(
        "--min_value",
        type=float,
        default=0,
        help="Minimum value of the elements initially",
    )
    parser.add_argument(
        "--max_range",
        type=float,
        default=np.inf,
        help="Maximum range of the elements initially",
    )
    parser.add_argument(
        "--min_range",
        type=float,
        default=-np.inf,
        help="Minimum range of the elements initially",
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
        "--probability", type=float, default=0.001, help="Probability of new generation"
    )
    parser.add_argument(
        "--func",
        type=str,
        default="sphere",
        choices=[
            "ackley",
            "beale",
            "booth",
            "bukin6",
            "crossintray",
            "easom",
            "eggholder",
            "goldstein",
            "himmelblau",
            "holdertable",
            "levi",
            "matyas",
            "rastringin",
            "rosenbrock",
            "schaffer2",
            "sphere",
            "threehump",
        ],
    )
    return parser.parse_args()


def get_distance(x, y):
    return np.max(np.sum((x - y) ** 2, axis=1)) ** 0.5


def main():
    # initialisation
    args = _parse_args()
    num_units = args.num_units
    num_feat = args.num_feat
    if args.func == "ackley":
        function = func.ackley
        args.max_range, args.min_range = (32, -32)
    elif args.func == "beale":
        function = func.beale
        args.max_range, args.min_range = (4.5, -4.5)
    elif args.func == "booth":
        function = func.booth
        args.max_range, args.min_range = (10, -10)
    elif args.func == "bukin6":
        function = func.bukin6
        args.max_range, args.min_range = (3, -3)
    elif args.func == "crossintray":
        function = func.crossintray
        args.max_range, args.min_range = (10, -10)
    elif args.func == "easom":
        function = func.easom
        args.max_range, args.min_range = (100, -100)
    elif args.func == "eggholder":
        function = func.eggholder
        args.max_range, args.min_range = (512, -512)
    elif args.func == "goldstein":
        function = func.goldstein
        args.max_range, args.min_range = (2, -2)
    elif args.func == "himmelblau":
        function = func.himmelblau
        args.max_range, args.min_range = (5, -5)
    elif args.func == "holdertable":
        function = func.holdertable
        args.max_range, args.min_range = (10, -10)
    elif args.func == "levi":
        function = func.levi
        args.max_range, args.min_range = (10, -10)
    elif args.func == "matyas":
        function = func.matyas
        args.max_range, args.min_range = (10, -10)
    elif args.func == "rastringin":
        function = func.rastrigin
        args.max_range, args.min_range = (5.12, -5.12)
    elif args.func == "rosenbrock":
        function = func.rosenbrock
    elif args.func == "schaffer2":
        function = func.schaffer2
        args.max_range, args.min_range = (100, -100)
    elif args.func == "sphere":
        function = func.sphere
    elif args.func == "threehump":
        function = func.threehump
        args.max_range, args.min_range = (5, -5)
    else:
        function = f
    position = np.random.random((num_units, num_feat))
    position = position * (args.max_value - args.min_value) + args.min_value
    best_individual = np.copy(position)
    best_pos = position[np.argmin(f(best_individual)), :].reshape(1, num_feat)
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
        if args.modified:
            prob = np.random.random((args.num_units, 1)) < args.probability
            position = position * (1 - prob) + prob * (
                best_pos
                + np.random.randn(*best_pos.shape)
                * get_distance(best_pos, position)
                / 100
            )
            velocity = velocity * (1 - prob)
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


if __name__ == "__main__":
    main()
