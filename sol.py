import numpy as np

from scipy.integrate import odeint
from Params import Params


class Solution:

    __end_time = None
    __start_time = 0
    __continuous_period = 0.001
    __KC_Continuous = None
    __KC_Discrete = None
    __exp_param = None
    __param = None

    def __init__(self, end_time, params: Params):
        self.__end_time = end_time
        self.__exp_param = -0.01
        self.__param = params.get_rho0()
        self.__KC_Continuous = params.calculate_KC_Continuous()
        self.__KC_Discrete = params.calculate_KC_Discrete()

    def influence_harmonic(self, t):
        return np.array([
            [np.sin(self.__param * t)],
            [np.sin(self.__param * t)],
            [np.sin(self.__param * t)],
            [np.sin(self.__param * t)],
            [np.sin(self.__param * t)],
            [np.sin(self.__param * t)]
        ])

    def influence_dumped(self, t):
        return 0.1 * np.exp(self.__exp_param * t) * np.array([
            [0],
            [0],
            [np.sin(t)],
            [0],
            [0],
            [np.cos(t)]
        ])

    def get_continuous_solution(self, params: Params, is_harmonic: bool):
        x0 = np.zeros([12, 1])
        time = np.arange(self.__start_time, self.__end_time + self.__continuous_period, self.__continuous_period)
        count = time.size
        z1_cont = np.zeros([2, count])
        z2_cont = np.zeros([4, count])
        z1_cont_not_norm = np.zeros([4, count])

        A = params.get_A()
        Bu = params.get_Bu()
        Bw = params.get_Bw()
        CL = params.get_CL()

        z2_cont[:, 0:1] = self.__KC_Continuous @ x0
        coord = CL @ x0
        z1_cont_not_norm[:, 0:1] = coord
        z1_cont[:, 0:1] = np.array([
            [np.linalg.norm(coord[0:3:2])],
            [np.linalg.norm(coord[1:4:2])]
        ])

        def vdp1(y, t, Ac, Bw):
            dydt = Ac @ y + (Bw @ self.influence_dumped(t)).reshape([12])
            return [dydt[0], dydt[1], dydt[2], dydt[3], dydt[4], dydt[5], dydt[6], dydt[7], dydt[8], dydt[9],
                                    dydt[10], dydt[11]]

        def vdp2(y, t, Ac, Bw):
            dydt = Ac @ y + (Bw @ self.influence_harmonic(t)).reshape([12])
            return [dydt[0], dydt[1], dydt[2], dydt[3], dydt[4], dydt[5], dydt[6], dydt[7], dydt[8], dydt[9],
                    dydt[10], dydt[11]]

        if is_harmonic:
            y = odeint(vdp2, x0.reshape([12]), time, args=(A + Bu @ self.__KC_Continuous, Bw))
        else:
            y = odeint(vdp1, x0.reshape([12]), time, args=(A + Bu @ self.__KC_Continuous, Bw))

        for k in range(count - 1):
            x0_ = np.transpose(y[k+1:k+2])
            z2_cont[:, k + 1:k + 2] = self.__KC_Continuous @ x0_
            coord = CL @ x0_
            z1_cont[:, k + 1:k + 2] = np.array([
                [np.linalg.norm(coord[0:3:2])],
                [np.linalg.norm(coord[1:4:2])]
            ])
            z1_cont_not_norm[:, k+1:k+2] = coord

        return z1_cont, z2_cont, z1_cont_not_norm

    def get_discrete_solution(self, params: Params, is_harmonic: bool):
        x0 = np.zeros([12, 1])
        step = params.get_step()
        time = np.arange(self.__start_time, self.__end_time + step, step)
        count = time.size
        z1_disc = np.zeros([2, count])
        z2_disc = np.zeros([4, count])
        z1_disc_not_norm = np.zeros([4, count])

        Ad = params.get_Ad()
        Bdu = params.get_Bdu()
        Bdw = params.get_Bdw()
        CL = params.get_CL()

        for k in range(count):
            z2_disc[:, k:k+1] = self.__KC_Discrete @ x0
            coord = CL @ x0
            z1_disc[:, k:k+1] = [
                [np.linalg.norm(coord[0:3:2])],
                [np.linalg.norm(coord[1:4:2])]
            ]
            z1_disc_not_norm[:, k:k+1] = coord
            if is_harmonic:
                x0 = (Ad + Bdu @ self.__KC_Discrete) @ x0 + Bdw @ self.influence_harmonic(time[k])
            else:
                x0 = (Ad + Bdu @ self.__KC_Discrete) @ x0 + Bdw @ self.influence_dumped(time[k])

        return z1_disc, z2_disc, z1_disc_not_norm

    def get_solution(self, params: Params, is_harmonic: bool):
        step = params.get_step()
        time = np.arange(self.__start_time, self.__end_time + self.__continuous_period, self.__continuous_period)
        count = time.size
        z1 = np.zeros([2, count])
        z2 = np.zeros([4, count])
        z1_not_norm = np.zeros([4, count])
        x0 = np.zeros([12, 1])

        A = params.get_A()
        Bu = params.get_Bu()
        Bw = params.get_Bw()
        CL = params.get_CL()

        z2[:, 0: 1] = self.__KC_Discrete @ x0
        coord = CL @ x0

        z1[:, 0:1] = np.array([
             [np.linalg.norm(coord[0:3:2])],
             [np.linalg.norm(coord[1:4:2])]
        ])
        z1_not_norm[:, 0:1] = coord

        def vdp1(y, t, _A, _Bwi, _BwKCx0):
            dydt = A @ y + _Bwi.reshape([12]) + _BwKCx0.reshape([12])
            return [dydt[0], dydt[1], dydt[2], dydt[3], dydt[4], dydt[5], dydt[6], dydt[7], dydt[8], dydt[9],
                    dydt[10], dydt[11]]

        for k in range(count - 1):
            if is_harmonic:
                y = odeint(vdp1, x0.reshape([12]), [time[k], time[k + 1]],
                           args=(A, Bw @ self.influence_harmonic(time[k + 1]),
                           Bu @ self.__KC_Discrete @ x0.reshape([12, 1])))
            else:
                y = odeint(vdp1, x0.reshape([12]), [time[k], time[k + 1]],
                           args=(A, Bw @ self.influence_dumped(time[k + 1]),
                           Bu @ self.__KC_Discrete @ x0.reshape([12, 1])))

            x0 = np.transpose(y[1:2])

            z2[:, k + 1:k + 2] = self.__KC_Discrete @ x0
            coord = CL @ x0
            z1[:, k + 1:k + 2] = np.array([
                 [np.linalg.norm(coord[0:3:2])],
                 [np.linalg.norm(coord[1:4:2])]
            ])
            z1_not_norm[:, k + 1:k + 2] = coord

        return z1, z2, z1_not_norm




