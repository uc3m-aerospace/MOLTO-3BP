import json

from Common.Controllers.OrbitController import OrbitController
import Utility.helpers as h
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import solve_ivp
from scipy import linalg
import os
import Utility.exceptions as exe
import Common.Config.configurator as conf

class LyapunovController(OrbitController):
    LYAPUNOV_ORBIT_SAMPLE_PATH = '\Externals\lyapunov_orbit_sample.txt'

    def __init__(self, data):
        super().__init__(attributes={"data": json.dumps(data)})
        if not self._data:
            self._data = data


    def main(self):
        '''
        Main implementation class for Input Manifold processor for Lyapunov
        :return:
        '''
        pass

    def generator(self):
        '''
        Main implementation class for computation
        :return:
        '''
        pass

    def plot(self):
        '''
        Plotting specificity implementation class
        :return:
        '''
        pass