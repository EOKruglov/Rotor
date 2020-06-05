import cvxpy as cvx
import math
import numpy as np
import scipy
import sympy
from scipy import integrate
from sympy import Matrix, symbols


class Params:

    _J = None
    _Jz = None
    _l0 = None
    _m = None
    _c = None
    _s0 = None
    _omega = None
    _T = None
    _rho0 = None
    _mu0 = None
    _lambda = None
    _eps = None
    _nu1 = None
    _nu2 = None
    _step = None
    _alpha = None

    _A = None
    _Bw = None
    _Bu = None

    _Ad = None
    _Bdu = None
    _Bdw = None

    _CL = None
    _DL = None
    _CR = None
    _DR = None

    _KC_Continious = None
    _KC_Discrete = None

    _B2 = None

    def __init__(self):
        pass

    def _set_J(self, J):
        self._J = J

    def _get_J(self):
        return self._J

    def _set_Jz(self, Jz):
        self._Jz = Jz

    def _get_Jz(self):
        return self._Jz

    def _set_l0(self, l0):
        self._l0 = l0

    def _get_l0(self):
        return self._l0

    def _set_m(self, m):
        self._m = m

    def _get_m(self):
        return self._m

    def _set_c(self, c):
        self._c = c

    def _get_c(self):
        return self._c

    def _set_s0(self, s0):
        self._s0 = s0

    def _get_s0(self):
        return self._s0

    def _set_omega(self, omega):
        self._omega = omega

    def _get_omega(self):
        return self._omega

    def _set_T(self, T):
        self._T = T

    def _get_T(self):
        return self._T

    def _set_rho0(self, rho0):
        self._rho0 = rho0

    def get_rho0(self):
        return self._rho0

    def _set_mu0(self, mu0):
        self._mu0 = mu0

    def _get_mu0(self):
        return self._mu0

    def _set_lambda(self, lambda0):
        self._lambda = lambda0

    def _get_lambda(self):
        return self._lambda

    def _set_eps(self, eps):
        self._eps = eps

    def _get_eps(self):
        return self._eps

    def _set_nu1(self, nu1):
        self._nu1 = nu1

    def _get_nu1(self):
        return self._nu1

    def _set_nu2(self, nu2):
        self._nu2 = nu2

    def _get_nu2(self):
        return self._nu2

    def _set_step(self, step):
        self._step = step

    def _get_step(self):
        return self._step

    def _set_alpha(self, alpha):
        self._alpha = alpha

    def _get_alpha(self):
        return self._alpha

    def _set_A(self, A):
        self._A = A

    def get_A(self):
        return self._A

    def _set_Bw(self, Bw):
        self._Bw = Bw

    def get_Bw(self):
        return self._Bw

    def _set_Bu(self, Bu):
        self._Bu = Bu

    def get_Bu(self):
        return self._Bu

    def _set_Ad(self, Ad):
        self._Ad = Ad

    def get_Ad(self):
        return self._Ad

    def _set_Bdu(self, Bdu):
        self._Bdu = Bdu

    def get_Bdu(self):
        return self._Bdu

    def _set_Bdw(self, Bdw):
        self._Bdw = Bdw

    def get_Bdw(self):
        return self._Bdw

    def _set_CL(self, CL):
        self._CL = CL

    def get_CL(self):
        return self._CL

    def _set_DL(self, DL):
        self._DL = DL

    def get_DL(self):
        return self._DL

    def _set_CR(self, CR):
        self._CR = CR

    def get_CR(self):
        return self._CR

    def _set_DR(self, DR):
        self._DR = DR

    def get_DR(self):
        return self._DR

    def _set_KC_Continious(self, KC_Continious):
        self._KC_Continious = KC_Continious

    def get_KC_Continious(self):
        return self._KC_Continious

    def _set_KC_Discrete(self, KC_Discrete):
        self._KC_Discrete = KC_Discrete

    def get_KC_Discrete(self):
        return self._KC_Discrete

    def _set_B2(self, B2):
        self._B2 = B2

    def _get_B2(self):
        return self._B2

    def initialize_input_variables(self, J, Jz, l0, m, c, s0, omega, eps, nu1, nu2, step, alpha):
        self._set_J(J)
        self._set_Jz(Jz)
        self._set_l0(l0)
        self._set_m(m)
        self._set_c(c)
        self._set_s0(s0)
        self._set_omega(omega)
        self._set_eps(eps)
        self._set_nu1(nu1)
        self._set_nu2(nu2)
        self._set_step(step)
        self._set_alpha(alpha)

    def initialize_non_input_variables(self):
        J = self._get_J()
        Jz = self._get_Jz()
        l0 = self._get_l0()
        m = self._get_m()
        c = self._get_c()
        s0 = self._get_s0()
        omega = self._get_omega()

        T = math.sqrt(m * s0)
        self._set_T(T)

        rho0 = (Jz * T * omega)/(J + m * (l0 ** 2))
        self._set_rho0(rho0)

        mu0 = (4 * m * (l0 ** 2)) / (J + m * (l0 ** 2))
        self._set_mu0(mu0)

        lambda0 = (c * (T ** 2)) / (J + m * (l0 ** 2))
        self._set_lambda(lambda0)

    def initialize_matrices(self):
        mu0 = self._get_mu0()
        eps = self._get_eps()
        rho0 = self.get_rho0()
        lambda0 = self._get_lambda()
        step = self._get_step()
        omega = self._get_omega()

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
            [-lambda0, lambda0, 0, -eps*omega, eps*omega, 0],
            [lambda0, -lambda0, 0, eps*omega, -eps*omega, 0],
            [0, 0, 0, 0, 0, 0],
            [eps*omega, -eps*omega, 0, -lambda0, lambda0, 0],
            [-eps*omega, eps*omega, 0, lambda0, -lambda0, 0],
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

        A = np.concatenate([
            np.concatenate([np.zeros([6, 6]), np.eye(6)], axis=1),
            np.concatenate([np.linalg.inv(A0) @ A2, np.linalg.inv(A0) @ A1], axis=1)
        ], axis=0)
        Bw = np.concatenate([np.zeros([6, 6]), np.linalg.inv(A0) @ B1], axis=0)
        Bu = np.concatenate([np.zeros([6, 4]), np.linalg.inv(A0) @ B2], axis=0)

        self._set_B2(B2)
        self._set_A(A)
        self._set_Bw(Bw)
        self._set_Bu(Bu)

        Ad = np.eye(12) + step * A
        Bdw = step * Bw
        Bdu = step * Bu

        self._set_Ad(Ad)
        self._set_Bdu(Bdu)
        self._set_Bdw(Bdw)

        CL = np.array([
            [0, 0, 0, 0.5, 0, 1, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0.5, 1, 0, 0, 0, 0, 0, 0],
            [-0.5, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, -0.5, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        ])
        DL = np.zeros([4, 4])
        CR = np.zeros([4, 12])
        DR = np.eye(4)

        self._set_CL(CL)
        self._set_DL(DL)
        self._set_CR(CR)
        self._set_DR(DR)

    def calculate_KC_Discrete(self):
        A = self.get_A()
        B2 = self._get_B2()
        Ad = self.get_Ad()
        Bdw = self.get_Bdw()
        Bdu = self.get_Bdu()
        CL = self.get_CL()
        DL = self.get_DL()
        CR = self.get_CR()
        DR = self.get_DR()
        alpha = self._get_alpha()

        nx = A.shape[0]
        nu = B2.shape[1]
        nzl = CL.shape[0]
        nzr = CR.shape[0]

        Y = cvx.Variable((nx, nx), symmetric=True)
        Z = cvx.Variable((nu, nx))
        gamma = cvx.Variable()

        first_constraint_matrix = cvx.bmat([
            [Y, Y @ Ad.T + Z.T @ Bdu.T, np.zeros([12, 6])],
            [Ad @ Y + Bdu @ Z, Y, Bdw],
            [np.zeros([6, 12]), Bdw.T, np.eye(6)]
        ])
        second_constraint_matrix = cvx.bmat([
            [Y, Y @ CL[0:3:2].T + Z.T @ DL[0:3:2].T],
            [CL[0:3:2] @ Y + DL[0:3:2] @ Z, ((alpha ** 2) * gamma) * np.eye(2)]
        ])
        third_constraint_matrix = cvx.bmat([
            [Y, Y @ CL[1:4:2].T + Z.T @ DL[1:4:2].T],
            [CL[1:4:2] @ Y + DL[1:4:2] @ Z, ((alpha ** 2) * gamma) * np.eye(2)]
        ])
        fourth_constraint_matrix = cvx.bmat([
            [Y, Y @ CR[0:1].T + Z.T @ DR[0:1].T],
            [CR[0:1] @ Y + DR[0:1] @ Z, (((1 - alpha) ** 2) * gamma) * np.eye(1)]
        ])
        fifth_constraint_matrix = cvx.bmat([
            [Y, Y @ CR[1:2].T + Z.T @ DR[1:2].T],
            [CR[1:2] @ Y + DR[1:2] @ Z, (((1 - alpha) ** 2) * gamma) * np.eye(1)]
        ])
        sixth_constraint_matrix = cvx.bmat([
            [Y, Y @ CR[2:3].T + Z.T @ DR[2:3].T],
            [CR[2:3] @ Y + DR[2:3] @ Z, (((1 - alpha) ** 2) * gamma) * np.eye(1)]
        ])
        seventh_constraint_matrix = cvx.bmat([
            [Y, Y @ CR[3:4].T + Z.T @ DR[3:4].T],
            [CR[3:4] @ Y + DR[3:4] @ Z, (((1 - alpha) ** 2) * gamma) * np.eye(1)]
        ])

        constraints = [first_constraint_matrix >> 0,
                       second_constraint_matrix >> 0,
                       third_constraint_matrix >> 0,
                       fourth_constraint_matrix >> 0,
                       fifth_constraint_matrix >> 0,
                       sixth_constraint_matrix >> 0,
                       seventh_constraint_matrix >> 0,
                       Y >> 0,
                       gamma >= 0]

        obj = cvx.Minimize(gamma)
        problem = cvx.Problem(obj, constraints)
        opt_val = problem.solve(verbose=True, solver='SCS')
        KC_Discrete = Z.value @ np.linalg.inv(Y.value)
        self._set_KC_Discrete(KC_Discrete)


    def caclucate_KC_Continious(self):
        A = self.get_A()
        Bw = self.get_Bw()
        Bu = self.get_Bu()
        CL = self.get_CL()
        DL = self.get_DL()
        CR = self.get_CR()
        DR = self.get_DR()
        alpha = self._get_alpha()

        nx = A.shape[0]
        nu = Bu.shape[1]
        nzl = CL.shape[0]
        nzr = CR.shape[0]

        Y = cvx.Variable((nx, nx), symmetric=True)
        Z = cvx.Variable((nu, nx))
        gamma = cvx.Variable()

        first_constraint_matrix = cvx.bmat([
            [Y @ A.T + A @ Y + Bu @ Z + Z.T @ Bu.T + Bw @ Bw.T]
        ])
        second_constraint_matrix = cvx.bmat([
            [Y, Y @ CL[0:3:2].T + Z.T @ DL[0:3:2].T],
            [CL[0:3:2] @ Y + DL[0:3:2] @ Z, alpha ** 2 * gamma * np.eye(2)]
        ])
        third_constraint_matrix = cvx.bmat([
            [Y, Y @ CL[1:4:2].T + Z.T @ DL[1:4:2].T],
            [CL[1:4:2] @ Y + DL[1:4:2] @ Z, alpha ** 2 * gamma * np.eye(2)]
        ])
        fourth_constraint_matrix = cvx.bmat([
            [Y, Y @ CR[0:1].T + Z.T @ DR[0:1].T],
            [CR[0:1] @ Y + DR[0:1] @ Z, (1 - alpha) ** 2 * gamma * np.eye(1)]
        ])
        fifth_constraint_matrix = cvx.bmat([
            [Y, Y @ CR[1:2].T + Z.T @ DR[1:2].T],
            [CR[1:2] @ Y + DR[1:2] @ Z, (1 - alpha) ** 2 * gamma * np.eye(1)]
        ])
        sixth_constraint_matrix = cvx.bmat([
            [Y, Y @ CR[2:3].T + Z.T @ DR[2:3].T],
            [CR[2:3] @ Y + DR[2:3] @ Z, (1 - alpha) ** 2 * gamma * np.eye(1)]
        ])
        seventh_constraint_matrix = cvx.bmat([
            [Y, Y @ CR[3:4].T + Z.T @ DR[3:4].T],
            [CR[3:4] @ Y + DR[3:4] @ Z, (1 - alpha) ** 2 * gamma * np.eye(1)]
        ])

        constraints = [first_constraint_matrix << 0,
                       second_constraint_matrix >> 0,
                       third_constraint_matrix >> 0,
                       fourth_constraint_matrix >> 0,
                       fifth_constraint_matrix >> 0,
                       sixth_constraint_matrix >> 0,
                       seventh_constraint_matrix >> 0,
                       Y >> 0,
                       gamma >= 0]

        obj = cvx.Minimize(gamma)
        problem = cvx.Problem(obj, constraints)
        opt_val = problem.solve(verbose=True, solver='SCS')
        KC_Continous = Z.value @ np.linalg.inv(Y.value)
        self._set_KC_Continious(KC_Continous)
