import numpy as np

from scipy.integrate import odeint
from utils import Params


class Solution:
    _x = None
    _z1 = None
    _z2 = None
    _z1_continious = None
    _z2_continious = None
    _z1_discrete = None
    _z2_discrete = None
    _start_time = None
    _end_time = None
    _period = None

    def __init__(self, start_time, end_time, period):
        self._start_time = start_time
        self._end_time = end_time
        self._period = period

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

    def _set_z1_continious(self, z1_continious):
        self._z1_continious = z1_continious

    def get_z1_continious(self):
        return self._z1_continious

    def _set_z2_continious(self, z2_continious):
        self._z2_continious = z2_continious

    def get_z2_continious(self):
        return self._z2_continious

    def _set_z1_discrete(self, z1_discrete):
        self._z1_discrete = z1_discrete

    def get_z1_discrete(self):
        return self._z1_discrete

    def _set_z2_discrete(self, z2_discrete):
        self._z2_discrete = z2_discrete

    def get_z2_discrete(self):
        return self._z2_discrete

    def influence_harmonic(self, t, rho0):
        return np.array([
            [np.sin(rho0 * t)],
            [np.sin(rho0 * t)],
            [np.sin(rho0 * t)],
            [np.sin(rho0 * t)],
            [np.sin(rho0 * t)],
            [np.sin(rho0 * t)]
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

    def get_solution(self, params: Params, is_harmonic: bool = False):
        KC_Discrete = params.get_KC_Discrete()
        A = params.get_A()
        Bw = params.get_Bw()
        Bu = params.get_Bu()
        time = np.arange(self._get_start_time(), self._get_end_time() + self._get_period(), self._get_period())
        x0 = np.zeros([12, 1])
        CL = params.get_CL()
        count = time.size
        z1 = np.zeros([2, count])
        z2 = np.zeros([4, count])
        rho0 = params.get_rho0()

        z2[:, 0: 1] = KC_Discrete @ x0
        coord = CL @ x0

        z1[:, 0:1] = np.array([
            [np.linalg.norm(coord[0:3:2])],
            [np.linalg.norm(coord[1:4:2])]
        ])

        for k in range(count - 1):

            def vdp1(y, t, _A, _Bwi, _BwKCx0):
                dydt = A @ y + _Bwi.reshape([12]) + _BwKCx0.reshape([12])
                return [dydt[0], dydt[1], dydt[2], dydt[3], dydt[4], dydt[5], dydt[6], dydt[7], dydt[8], dydt[9],
                        dydt[10], dydt[11]]

            if is_harmonic:
                y = odeint(vdp1, x0.reshape([12]), [time[k], time[k + 1]], args=(A, Bw @ self.influence_harmonic(time[k], rho0),
                                                                                 Bu @ KC_Discrete @ x0.reshape(
                                                                                     [12, 1])))
            else:
                y = odeint(vdp1, x0.reshape([12]), [time[k], time[k + 1]], args=(A, Bw @ self.influence_dumped(time[k]),
                                                                                 Bu @ KC_Discrete @ x0.reshape([12, 1])))
            x0 = np.transpose(y[1:2])

            z2[:, k + 1:k + 2] = KC_Discrete @ x0
            coord = CL @ x0
            z1[:, k + 1:k + 2] = np.array([
                [np.linalg.norm(coord[0:3:2])],
                [np.linalg.norm(coord[1:4:2])]
            ])

        self._set_z1(z1)
        self._set_z2(z2)

        x0 = np.zeros([12, 1])
        KC_Continious = params.get_KC_Continious()
        z1_cont = np.zeros([2, count])
        z2_cont = np.zeros([4, count])

        z2_cont[:, 0:1] = KC_Continious @ x0
        coord = CL @ x0
        z1_cont[:, 0:1] = np.array([
            [np.linalg.norm(coord[0:3:2])],
            [np.linalg.norm(coord[1:4:2])]
        ])

        def vdp3(y, t, Ac, Bw):
            dydt = Ac @ y + (Bw @ self.influence_dumped(t)).reshape([12])
            return [dydt[0], dydt[1], dydt[2], dydt[3], dydt[4], dydt[5], dydt[6], dydt[7], dydt[8], dydt[9],
                                    dydt[10], dydt[11]]

        def vdp4(y, t, Ac, Bw, _rho0):
            dydt = Ac @ y + (Bw @ self.influence_harmonic(t, _rho0)).reshape([12])
            return [dydt[0], dydt[1], dydt[2], dydt[3], dydt[4], dydt[5], dydt[6], dydt[7], dydt[8], dydt[9],
                    dydt[10], dydt[11]]

        if is_harmonic:
            y = odeint(vdp4, x0.reshape([12]), time, args=(A + Bu @ KC_Continious, Bw, rho0))
        else:
            y = odeint(vdp3, x0.reshape([12]), time, args=(A + Bu @ KC_Continious, Bw))

        for k in range(count - 1):
            x0_ = np.transpose(y[k+1:k+2])
            z2_cont[:, k + 1:k + 2] = KC_Continious @ x0_
            coord = CL @ x0_
            z1_cont[:, k + 1:k + 2] = np.array([
                [np.linalg.norm(coord[0:3:2])],
                [np.linalg.norm(coord[1:4:2])]
            ])

        self._set_z1_continious(z1_cont)
        self._set_z2_continious(z2_cont)

        x0 = np.zeros([12, 1])
        z1_disc = np.zeros([2, count])
        z2_disc = np.zeros([4, count])
        Ad = params.get_Ad()
        Bdu = params.get_Bdu()
        Bdw = params.get_Bdw()

        for k in range(count):
            z2_disc[:, k:k+1] = KC_Discrete @ x0
            coord = CL @ x0
            z1_disc[:, k:k+1] = [
                [np.linalg.norm(coord[0:3:2])],
                [np.linalg.norm(coord[1:4:2])]
            ]
            if is_harmonic:
                x0 = (Ad + Bdu @ KC_Discrete) @ x0 + Bdw @ self.influence_harmonic(time[k], rho0)
            else:
                x0 = (Ad + Bdu @ KC_Discrete) @ x0 + Bdw @ self.influence_dumped(time[k])

        self._set_z1_discrete(z1_disc)
        self._set_z2_discrete(z2_disc)
