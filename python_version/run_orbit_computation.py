import argparse
import json
import math
import multiprocessing
import random

from decorator import contextmanager

import Utility.helpers as h

import numpy as np
import Common.Config.configurator as c
from Common.Controllers.ManifoldController import ManifoldController
from Common.Service.benchmark import Benchmark as b

class  TrajectoryCompute:
    '''
    Main Functionalities (expected string)
    Type of orbit of interest,
    'LY' == lyapunov orbit (2D)
    'HL' == halo orbit (3D)
    defaults to random
    '''

    ## type of Orbit to compute
    # Defaults to random choice
    type = random.choice(['HL','LY'])

    # System analyzed
    # Accepted values: 'SE', 'EM'
    # Sun-Earth, Earth-Moon, Whatever-CelestialBody, Saturn-(1 of its Moons)
    f = 'EM'

    # Case choice (Orbital amplitude)
    # Lyapunov Orbit characterization
    # Only x needed (0 < x < 1e-2)
    Ax = 1.05e-3

    # Halo orbit amplitude (perpendicular to the plane)
    # Halo Orbit characterization (a big number because its not nondimensional/not normalised)
    # In km, not normalised
    # TODO: maybe have a length agnostic normalised versioned
    Az = 110e3

    # 1 -> Northern variant of the orbit
    # Rotation constant
    # Easiest is no rotation [-PI/2; +PI/2]
    phi = 0

    # 1 or 3 -> Southern variant of the orbit
    m = 1

    # 1 -> Lagrange point L1
    # 2 -> Lagrange point L2
    # Orbit will sit around these 2 points [1,2] in |N
    LP = 2

    # Angle (in degrees) between the section required and
    # the +X semiplane taken from the 2nd primary
    poincSec = np.array([-90, 90])


    # Numerical parameters
    # Number of points in the orbit to propagate manifolds from
    # Linear progression on the 3D plot of blue lines
    # Suggested 3 < n < 9
    npoints = 5

    # The program propagates the perturbation in the direction chosen
    # both directions 0, interior realm 1, exterior realm -1
    d = 1

    # Propagation of stable branch, unstable one or both
    # both branches 0, unstable branch -1, stable branch 1
    # 0 => BLUE, RED (Towards, leaving)
    # 1 => BLUE (Towards the orbit)
    # -1 => RED (Leaving)
    # TODO Red or Blue pill in the docs
    branch = 0

    # print time period
    # discretisation of the trajectories
    # granularity/specificity (points in continuous non linear lines)
    prnt_out_dt = 0.001

    # Sets the program to evaluate a mission from orbit 1 .. n
    # Following intersections in 1 .. n Poincaré sections
    # The only values overwritten by this parameter are the case
    # choice variables and d and branch parameters to create the chain
    # Linear computational bloating
    # 0 for only one type, 1 for file reading
    input_seq = 0

    # Enables the introduction of data via text instead of typing
    text = 1

    # In order to introduce data from the exterior, the program expects a .txt input
    file = '/Externals/Samples/invariant_manifolds_sample_second.txt'

    parallel = 1

    def __init__(self, args={}):
        '''
        Initialize main features
        :param args:
        '''

        if len([v is not None for k, v in args.items()]) > 0:
            if 'f' in args and args.get("mode") is not None:
                self.f = args.get("mode")
            if 'Ax' in args and args.get("Ax") is not None:
                self.Ax = float(args.get("Ax"))
            if 'Az' in args and args.get("Az") is not None:
                self.Az = float(args.get("Az"))
            if 'phi' in args and args.get("phi") is not None:
                self.phi = float(args.get("phi"))
            if 'm' in args and args.get("m") is not None:
                self.m = int(args.get("m"))
            if 'LP' in args and args.get("LP") is not None:
                self.LP = int(args.get("LP"))
            if 'poincSec' in args and args.get("pointSec") is not None:
                self.poincSec = np.array(args.get("poincSec"))

            if 'n_points' in args and args.get("n_points") is not None:
                self.npoints = int(args.get("npoints"))

            if 'd' in args and args.get("d") is not None:
                self.d = int(args.get("d"))


            if 'branch' in args and args.get("branch") is not None:
                self.d = int(args.get("branch"))

            if 'prnt_out_dt' in args and args.get("print_out_dt") is not None:
                self.prnt_out_dt = float(args.get("prnt_out_dt"))

            if 'input_seq' in args and args.get("input_seq") is not None:
                self.prnt_out_dt = int(args.get("input_sq"))

            if 'text' in args and args.get("text") is not None:
                self.text = int(args.get("text"))

            if 'file' in args and args.get("file") is not None:
                self.file = args.get("file")

            if 'parallel' in args and args.get("parallel") is not None:
                self.parallel = args.get("parallel")


    @contextmanager
    def poolcontext(*args, **kwargs):
        """
        A bit of context never hurts....
        """

        pool = multiprocessing.Pool(*args, **kwargs)
        yield pool
        pool.terminate()

    def main(self):
        input = {'type': self.type, 'mode': self.f, 'Ax_tgt': self.Ax, 'Az': self.Az, 'm': self.m, 'phi': self.phi,
                 'LP': self.LP, 'poincSec': self.poincSec, 'npoints': self.npoints, 'd': self.d, 'branch': self.branch,
                 'prnt_out_dt': self.prnt_out_dt, 'input_seq': self.input_seq, 'text': self.text, 'file': self.file}

        if self.parallel > 1:
            # TODO change processes here for split plot or orbits
            config = c.Configurator()
            core_divisor = config.get('run.core_split', 3)
            max_async_process = math.ceil(multiprocessing.cpu_count() / core_divisor)
            with multiprocessing.Pool(processes=max_async_process) as pool:
                # Partial over agents (as we run all those selected accounts as per one optimisation)
                p = pool.map(ManifoldController.main_runner(data=input), self.parallel)
                # close the pool and wait for the work to finish
                pool.close()
                pool.join()
        else:
            ManifoldController.main_runner(data=input)





if __name__ == '__main__':


    parser = argparse.ArgumentParser(
        description='Sync of the accounts')
    parser.add_argument('-type', '--type', type=str, help='', required=False)
    parser.add_argument('-f', '--mode',  type=str, help='', required=False)
    parser.add_argument('-Ax', '--Ax', type=float, help='',  required=False)
    parser.add_argument('-Az', '--Az', type=float, help='',  required=False)
    parser.add_argument('-m', '--m', type=int, required=False)
    parser.add_argument('-phi', '--phi', type=float, required=False)
    parser.add_argument('-LP', '--LP', type=int, required=False)

    parser.add_argument('-poincSec', '--poincSec', type=int, help='',  required=False)

    parser.add_argument('-npoints', '--npoints',  type=float, help='',  required=False)

    parser.add_argument('-d', '--d',  type=int,
                        help='',  required=False)
    parser.add_argument('-branch', '--branch',  type=int,
                        help='',  required=False)
    parser.add_argument('-print_out_dt', '--print_out_dt',  type=float,
                        help='',  required=False)
    parser.add_argument('-input_seq', '--input_seq',  type=int,
                        help='',  required=False)
    parser.add_argument('-text', '--text',  type=int,
                        help='',  required=False)
    parser.add_argument('-file', '--file',  type=str,
                        help='',  required=False)

    args = parser.parse_args()
    h.log(f'----Started Orbit Trajectory compute process for {args}----')

    b.print_stack(f'----Started Orbit Trajectory compute process for {args}----')

    # TODO, perhaps unwrap here

    TrajectoryCompute(
        args=vars(args)
    ).main()

    b.print_stack(f'----Finished Orbit Trajectory compute with args {args}----')
