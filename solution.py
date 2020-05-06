import numpy as np

from scipy.integrate import odeint
from utils import Params


class Solution:
    _x = None
    _z1 = None
    _z2 = None
    _start_time = None
    _end_time = None
    _period = None

    def __init__(self, start_time, end_time, period):
        _start_time = start_time
        _end_time = end_time
        _period = period

    def _get_start_time(self):
        return self._start_time

    def _get_end_time(self):
        return self._end_time

    def _get_period(self):
        return self._period

    def _set_x(self, x):
        self._x = x

    def _get_x(self):
        return self._x

    def _set_z1(self, z1):
        self._z1 = z1

    def get_z1(self):
        return self._z1

    def _set_z2(self, z2):
        self._z2 = z2

    def get_z2(self):
        return self._z2

    def influence_harmonic(self, t, rho0):
        return np.array([
            [np.sin(rho0 * t)],
            [np.sin(rho0 * t)],
            [np.sin(rho0 * t)],
            [np.sin(rho0 * t)],
            [np.sin(rho0 * t)],
            [np.sin(rho0 * t)0]
        ])

    def influence_dumped(self, t):
        return 0.1 * np.exp(-0.01 * t) * np.array([
            [0],
            [0],
            [np.sin(t)],
            [0],
            [0],
            [np.cos(t)]
        ])


    def get_solution(self, params: Params):
        KC_Discrete = params.get_KC_Discrete()
        A = params.get_A()
        Bw = params.get_Bw()
        Bu = params.get_Bu()
        t = np.arange(self._get_start_time(), self._get_end_time() + self._get_period(), self._get_period())
        x0 = np.zeros([12, 1])
        CL = params.get_CL()
        count = t.size
        z1 = np.zeros([2, count])
        z2 = np.zeros([4, count])

        z2[:, 0] = KC_Discrete * x0
        coord = CL * x0

        z1[:, 0] = np.array([
            [np.linalg.norm(coord[0:3:2])],
            [np.linalg.norm(coord[1:4:2])]
        ])

        for k in range(count - 1):
            xg = odeint(A * xg + Bw * self.influence_dumped(t[k]) + Bu * KC_Discrete * x0, x0,
                        np.arange(t[k], t[k + 1], self._get_period()))

            x0 = np.transpose(xg[-1, :])

            z2[:, k + 1] = KC_Discrete * x0
            coord = CL * x0
            z1[:, k + 1] = np.array([
            [np.linalg.norm(coord[0:3:2])],
            [np.linalg.norm(coord[1:4:2])]
        ])

        self._set_z1(z1)
        self._set_z2(z2)

