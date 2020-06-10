import cvxpy as cvx
import numpy as np


class Params:

    __J = None
    __Jz = None
    __l0 = None
    __m = None
    __c = None
    __s0 = None
    __omega = None
    __eps = None
    __step = None
    __alpha = None

    __A = None
    __Bw = None
    __Bu = None

    __Ad = None
    __Bdu = None
    __Bdw = None

    __CL = None
    __DL = None
    __CR = None
    __DR = None

    __B2 = None
    __rho0 = None

    def __init__(self, J, Jz, l0, m, c, s0, omega, eps, step, alpha):
        self.__J = J
        self.__Jz = Jz
        self.__l0 = l0
        self.__m = m
        self.__c = c
        self.__s0 = s0
        self.__omega = omega
        self.__eps = eps
        self.__step = step
        self.__alpha = alpha

        T = np.sqrt(m * s0)
        mu0 = (4 * m * (l0 ** 2)) / (J + m * (l0 ** 2))
        rho0 = (Jz * T * omega) / (J + m * (l0 ** 2))
        self.__rho0 = rho0
        lambda0 = (c * (T ** 2)) / (J + m * (l0 ** 2))

        A0 = np.array([
            [1, 0, -mu0, 0, 0, 0],
            [0, 1, mu0, 0, 0, 0],
            [-0.125, 0.125, 1, 0, 0, 0],
            [0, 0, 0, 1, 0, mu0],
            [0, 0, 0, 0, 1, -mu0],
            [0, 0, 0, 0.125, -0.125, 1]
        ])

        A1 = np.array([
            [-eps, eps, 0, -rho0, 0, 0],
            [eps, -eps, 0, 0, -rho0, 0],
            [0, 0, 0, 0, 0, 0],
            [rho0, 0, 0, -eps, eps, 0],
            [0, rho0, 0, eps, -eps, 0],
            [0, 0, 0, 0, 0, 0]
        ])

        A2 = np.array([
            [-lambda0, lambda0, 0, -eps * omega, eps * omega, 0],
            [lambda0, -lambda0, 0, eps * omega, -eps * omega, 0],
            [0, 0, 0, 0, 0, 0],
            [eps * omega, -eps * omega, 0, -lambda0, lambda0, 0],
            [-eps * omega, eps * omega, 0, lambda0, -lambda0, 0],
            [0, 0, 0, 0, 0, 0]
        ])

        B1 = np.eye(6)

        B2 = np.array([
            [-4 * mu0, 0, 0, 0],
            [0, 4 * mu0, 0, 0],
            [1, 1, 0, 0],
            [0, 0, 4 * mu0, 0],
            [0, 0, 0, -4 * mu0],
            [0, 0, 1, 1]
        ])

        self.__B2 = B2

        A = np.concatenate([
            np.concatenate([np.zeros([6, 6]), np.eye(6)], axis=1),
            np.concatenate([np.linalg.inv(A0) @ A2, np.linalg.inv(A0) @ A1], axis=1)
        ], axis=0)
        Bw = np.concatenate([np.zeros([6, 6]), np.linalg.inv(A0) @ B1], axis=0)
        Bu = np.concatenate([np.zeros([6, 4]), np.linalg.inv(A0) @ B2], axis=0)

        self.__A = A
        self.__Bw = Bw
        self.__Bu = Bu

        Ad = np.eye(12) + step * A
        Bdw = step * Bw
        Bdu = step * Bu

        self.__Ad = Ad
        self.__Bdw = Bdw
        self.__Bdu = Bdu

        CL = np.array([
            [0, 0, 0, 2 * l0, 0, 1, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, -2 * l0, 1, 0, 0, 0, 0, 0, 0],
            [-2 * l0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 2 * l0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        ])
        DL = np.zeros([4, 4])
        CR = np.zeros([4, 12])
        DR = np.eye(4)

        self.__CL = CL
        self.__DL = DL
        self.__CR = CR
        self.__DR = DR

    def get_A(self):
        return self.__A

    def get_Bu(self):
        return self.__Bu

    def get_Bw(self):
        return self.__Bw

    def get_Ad(self):
        return self.__Ad

    def get_Bdu(self):
        return self.__Bdu

    def get_Bdw(self):
        return self.__Bdw

    def get_CL(self):
        return self.__CL

    def get_step(self):
        return self.__step

    def get_rho0(self):
        return self.__rho0

    def calculate_KC_Continuous(self):
        nx = self.__A.shape[0]
        nu = self.__Bu.shape[1]

        Y = cvx.Variable((nx, nx), symmetric=True)
        Z = cvx.Variable((nu, nx))
        gamma = cvx.Variable()

        first_constraint_matrix = cvx.bmat([
            [Y @ self.__A.T + self.__A @ Y + self.__Bu @ Z + Z.T @ self.__Bu.T, self.__Bw],
            [self.__Bw.T, -np.eye(6)]
        ])
        second_constraint_matrix = cvx.bmat([
            [Y, Y @ self.__CL[0:3:2].T + Z.T @ self.__DL[0:3:2].T],
            [self.__CL[0:3:2] @ Y + self.__DL[0:3:2] @ Z, self.__alpha ** 2 * gamma * np.eye(2)]
        ])
        third_constraint_matrix = cvx.bmat([
            [Y, Y @ self.__CL[1:4:2].T + Z.T @ self.__DL[1:4:2].T],
            [self.__CL[1:4:2] @ Y + self.__DL[1:4:2] @ Z, self.__alpha ** 2 * gamma * np.eye(2)]
        ])
        fourth_constraint_matrix = cvx.bmat([
            [Y, Y @ self.__CR[0:1].T + Z.T @ self.__DR[0:1].T],
            [self.__CR[0:1] @ Y + self.__DR[0:1] @ Z, (1 - self.__alpha) ** 2 * gamma * np.eye(1)]
        ])
        fifth_constraint_matrix = cvx.bmat([
            [Y, Y @ self.__CR[1:2].T + Z.T @ self.__DR[1:2].T],
            [self.__CR[1:2] @ Y + self.__DR[1:2] @ Z, (1 - self.__alpha) ** 2 * gamma * np.eye(1)]
        ])
        sixth_constraint_matrix = cvx.bmat([
            [Y, Y @ self.__CR[2:3].T + Z.T @ self.__DR[2:3].T],
            [self.__CR[2:3] @ Y + self.__DR[2:3] @ Z, (1 - self.__alpha) ** 2 * gamma * np.eye(1)]
        ])
        seventh_constraint_matrix = cvx.bmat([
            [Y, Y @ self.__CR[3:4].T + Z.T @ self.__DR[3:4].T],
            [self.__CR[3:4] @ Y + self.__DR[3:4] @ Z, (1 - self.__alpha) ** 2 * gamma * np.eye(1)]
        ])

        constraints = [first_constraint_matrix << 0,
                       second_constraint_matrix >> 0,
                       third_constraint_matrix >> 0,
                       fourth_constraint_matrix >> 0,
                       fifth_constraint_matrix >> 0,
                       sixth_constraint_matrix >> 0,
                       seventh_constraint_matrix >> 0,
                       Y >> 0,
                       gamma >= np.finfo(float).eps]

        obj = cvx.Minimize(gamma)
        problem = cvx.Problem(obj, constraints)
        opt_val = problem.solve(verbose=True, solver='SCS')
        KC_Continuous = Z.value @ np.linalg.inv(Y.value)
        return KC_Continuous


    def calculate_KC_Discrete(self):
        nx = self.__A.shape[0]
        nu = self.__B2.shape[1]

        Y = cvx.Variable((nx, nx), symmetric=True)
        Z = cvx.Variable((nu, nx))
        gamma = cvx.Variable()

        first_constraint_matrix = cvx.bmat([
            [Y, Y @ self.__Ad.T + Z.T @ self.__Bdu.T, np.zeros([12, 6])],
            [self.__Ad @ Y + self.__Bdu @ Z, Y, self.__Bdw],
            [np.zeros([6, 12]), self.__Bdw.T, np.eye(6)]
        ])
        second_constraint_matrix = cvx.bmat([
            [Y, Y @ self.__CL[0:3:2].T + Z.T @ self.__DL[0:3:2].T],
            [self.__CL[0:3:2] @ Y + self.__DL[0:3:2] @ Z, ((self.__alpha ** 2) * gamma) * np.eye(2)]
        ])
        third_constraint_matrix = cvx.bmat([
            [Y, Y @ self.__CL[1:4:2].T + Z.T @ self.__DL[1:4:2].T],
            [self.__CL[1:4:2] @ Y + self.__DL[1:4:2] @ Z, ((self.__alpha ** 2) * gamma) * np.eye(2)]
        ])
        fourth_constraint_matrix = cvx.bmat([
            [Y, Y @ self.__CR[0:1].T + Z.T @ self.__DR[0:1].T],
            [self.__CR[0:1] @ Y + self.__DR[0:1] @ Z, (((1 - self.__alpha) ** 2) * gamma) * np.eye(1)]
        ])
        fifth_constraint_matrix = cvx.bmat([
            [Y, Y @ self.__CR[1:2].T + Z.T @ self.__DR[1:2].T],
            [self.__CR[1:2] @ Y + self.__DR[1:2] @ Z, (((1 - self.__alpha) ** 2) * gamma) * np.eye(1)]
        ])
        sixth_constraint_matrix = cvx.bmat([
            [Y, Y @ self.__CR[2:3].T + Z.T @ self.__DR[2:3].T],
            [self.__CR[2:3] @ Y + self.__DR[2:3] @ Z, (((1 - self.__alpha) ** 2) * gamma) * np.eye(1)]
        ])
        seventh_constraint_matrix = cvx.bmat([
            [Y, Y @ self.__CR[3:4].T + Z.T @ self.__DR[3:4].T],
            [self.__CR[3:4] @ Y + self.__DR[3:4] @ Z, (((1 - self.__alpha) ** 2) * gamma) * np.eye(1)]
        ])

        constraints = [first_constraint_matrix >> 0,
                       second_constraint_matrix >> 0,
                       third_constraint_matrix >> 0,
                       fourth_constraint_matrix >> 0,
                       fifth_constraint_matrix >> 0,
                       sixth_constraint_matrix >> 0,
                       seventh_constraint_matrix >> 0,
                       Y >> 0,
                       gamma >= np.finfo(float).eps]

        obj = cvx.Minimize(gamma)
        problem = cvx.Problem(obj, constraints)
        opt_val = problem.solve(verbose=True, solver='SCS')
        KC_Discrete = Z.value @ np.linalg.inv(Y.value)
        return KC_Discrete
