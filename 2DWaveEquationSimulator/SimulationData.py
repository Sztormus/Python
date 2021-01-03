import math
import numpy as np
from scipy.linalg import lu_factor, lu_solve


class SimulationData:
    def __init__(self):
        self.__n_x = 25
        self.__n_y = 25
        self.__t = 1.0
        self.__alpha = 0.0
        self.__beta = 0.25
        self.__gamma = 0.5
        self.__fixed_boundary_conditions = True
        self.__zero_diaphragm_displacement = False
        self.__forcing_term_check = False

        self.__delta_x = 1.0 / (self.__n_x - 1)
        self.__delta_y = 1.0 / (self.__n_y - 1)
        self.__delta_t = 0.01
        self.__n_t = int(self.__t / self.__delta_t)
        self.__t_f = 0.5
        self.__x_f = self.__delta_x * int(self.__n_x / 2)
        self.__y_f = self.__delta_y * int(self.__n_y / 2)
        self.__sigma = 0.2

        self.__parameter_a = 0.5 - self.__gamma + self.__beta
        self.__parameter_b = 0.5 + self.__gamma - 2.0 * self.__beta
        self.__parameter_c = 1.0 + 2.0 * self.__delta_t ** 2.0 * self.__beta * (1.0 / self.__delta_x ** 2.0 + 1.0 / self.__delta_y ** 2.0 + self.__alpha / self.__delta_t)
        self.__parameter_d = self.__delta_t ** 2.0 * self.__beta / self.__delta_x ** 2.0
        self.__parameter_e = self.__delta_t ** 2.0 * self.__beta / self.__delta_y ** 2.0
        self.__parameter_f = 2.0 * (1.0 - self.__delta_t * self.__parameter_a * self.__alpha - self.__delta_t ** 2.0 * self.__parameter_b * (1.0 / self.__delta_x ** 2.0 + 1.0 / self.__delta_y ** 2.0 + self.__alpha / self.__delta_t) + self.__delta_t * self.__beta * self.__alpha)
        self.__parameter_g = 1.0 + 2.0 * self.__delta_t ** 2.0 * self.__parameter_a * (1.0 / self.__delta_x ** 2.0 + 1.0 / self.__delta_y ** 2.0 - self.__alpha / self.__delta_t) - 2.0 * self.__delta_t * self.__parameter_b * self.__alpha
        self.__parameter_h = self.__delta_t ** 2.0 * self.__parameter_a / self.__delta_x ** 2.0
        self.__parameter_j = self.__delta_t ** 2.0 * self.__parameter_a / self.__delta_y ** 2.0
        self.__parameter_k = self.__delta_t ** 2.0 * self.__parameter_b / self.__delta_x ** 2.0
        self.__parameter_l = self.__delta_t ** 2.0 * self.__parameter_b / self.__delta_y ** 2.0

        self.__s = np.zeros((self.__n_x * self.__n_y, self.__n_x * self.__n_y), dtype=np.float64)
        self.__b = np.zeros(self.__n_x * self.__n_y, dtype=np.float64)

        self.__u_previous = np.zeros(self.__n_x * self.__n_y, dtype=np.float64)
        self.__u_actual = np.zeros(self.__n_x * self.__n_y, dtype=np.float64)
        self.__u_next = np.zeros(self.__n_x * self.__n_y, dtype=np.float64)

        if not self.__zero_diaphragm_displacement:
            for k in range(self.__n_x * self.__n_y):
                j = math.floor(k / self.__n_x)
                i = k - j * self.__n_x
                if not (i == 0 or i == self.__n_x - 1 or j == 0 or j == self.__n_y - 1):
                    self.__u_previous[k] = (i * self.__delta_x - (i * self.__delta_x) ** 2) * (j * self.__delta_y - (j * self.__delta_y) ** 2)
                    self.__u_actual[k] = (i * self.__delta_x - (i * self.__delta_x) ** 2) * (j * self.__delta_y - (j * self.__delta_y) ** 2)

        if not self.__fixed_boundary_conditions:
            for k in range(self.__n_x * self.__n_y):
                j = math.floor(k / self.__n_x)
                i = k - j * self.__n_x
                if i == 0:
                    self.__s[k, k] = 1
                    self.__s[k, k + 1] = -1
                elif i == self.__n_x - 1:
                    self.__s[k, k] = 1
                    self.__s[k, k - 1] = -1
                elif j == 0:
                    self.__s[k, k] = 1
                    self.__s[k, k + self.__n_x] = -1
                elif j == self.__n_y - 1:
                    self.__s[k, k] = 1
                    self.__s[k, k - self.__n_x] = -1
                else:
                    self.__s[k, k] = self.__parameter_c
                    self.__s[k, k + 1] = self.__s[k, k - 1] = -self.__parameter_d
                    self.__s[k, k + self.__n_x] = self.__s[k, k - self.__n_x] = -self.__parameter_e
        else:
            for k in range(self.__n_x * self.__n_y):
                j = math.floor(k / self.__n_x)
                i = k - j * self.__n_x
                if i == 0 or i == self.__n_x - 1 or j == 0 or j == self.__n_y - 1:
                    self.__s[k, k] = 1
                else:
                    self.__s[k, k] = self.__parameter_c
                    self.__s[k, k + 1] = self.__s[k, k - 1] = -self.__parameter_d
                    self.__s[k, k + self.__n_x] = self.__s[k, k - self.__n_x] = -self.__parameter_e

    def change_initial_configuration(self, n_x, n_y, t, alpha, beta, gamma, fixed_boundary_conditions, zero_diaphragm_displacement, forcing_term_check):
        self.__n_x = n_x
        self.__n_y = n_y
        self.__t = t
        self.__alpha = alpha
        self.__beta = beta
        self.__gamma = gamma
        self.__fixed_boundary_conditions = fixed_boundary_conditions
        self.__zero_diaphragm_displacement = zero_diaphragm_displacement
        self.__forcing_term_check = forcing_term_check

        self.__delta_x = 1.0 / (self.__n_x - 1)
        self.__delta_y = 1.0 / (self.__n_y - 1)
        self.__delta_t = 0.01
        self.__n_t = int(self.__t / self.__delta_t)
        self.__t_f = 0.5
        self.__x_f = self.__delta_x * int(self.__n_x / 2)
        self.__y_f = self.__delta_y * int(self.__n_y / 2)
        self.__sigma = 0.2

        self.__parameter_a = 0.5 - self.__gamma + self.__beta
        self.__parameter_b = 0.5 + self.__gamma - 2.0 * self.__beta
        self.__parameter_c = 1.0 + 2.0 * self.__delta_t ** 2.0 * self.__beta * (1.0 / self.__delta_x ** 2.0 + 1.0 / self.__delta_y ** 2.0 + self.__alpha / self.__delta_t)
        self.__parameter_d = self.__delta_t ** 2.0 * self.__beta / self.__delta_x ** 2.0
        self.__parameter_e = self.__delta_t ** 2.0 * self.__beta / self.__delta_y ** 2.0
        self.__parameter_f = 2.0 * (1.0 - self.__delta_t * self.__parameter_a * self.__alpha - self.__delta_t ** 2.0 * self.__parameter_b * (1.0 / self.__delta_x ** 2.0 + 1.0 / self.__delta_y ** 2.0 + self.__alpha / self.__delta_t) + self.__delta_t * self.__beta * self.__alpha)
        self.__parameter_g = 1.0 + 2.0 * self.__delta_t ** 2.0 * self.__parameter_a * (1.0 / self.__delta_x ** 2.0 + 1.0 / self.__delta_y ** 2.0 - self.__alpha / self.__delta_t) - 2.0 * self.__delta_t * self.__parameter_b * self.__alpha
        self.__parameter_h = self.__delta_t ** 2.0 * self.__parameter_a / self.__delta_x ** 2.0
        self.__parameter_j = self.__delta_t ** 2.0 * self.__parameter_a / self.__delta_y ** 2.0
        self.__parameter_k = self.__delta_t ** 2.0 * self.__parameter_b / self.__delta_x ** 2.0
        self.__parameter_l = self.__delta_t ** 2.0 * self.__parameter_b / self.__delta_y ** 2.0

        self.__s = np.zeros((self.__n_x * self.__n_y, self.__n_x * self.__n_y), dtype=np.float64)
        self.__b = np.zeros(self.__n_x * self.__n_y, dtype=np.float64)

        self.__u_previous = np.zeros(self.__n_x * self.__n_y, dtype=np.float64)
        self.__u_actual = np.zeros(self.__n_x * self.__n_y, dtype=np.float64)
        self.__u_next = np.zeros(self.__n_x * self.__n_y, dtype=np.float64)

        if not self.__zero_diaphragm_displacement:
            for k in range(self.__n_x * self.__n_y):
                j = math.floor(k / self.__n_x)
                i = k - j * self.__n_x
                if not (i == 0 or i == self.__n_x - 1 or j == 0 or j == self.__n_y - 1):
                    self.__u_previous[k] = (i * self.__delta_x - (i * self.__delta_x) ** 2) * (j * self.__delta_y - (j * self.__delta_y) ** 2)
                    self.__u_actual[k] = (i * self.__delta_x - (i * self.__delta_x) ** 2) * (j * self.__delta_y - (j * self.__delta_y) ** 2)

        if not self.__fixed_boundary_conditions:
            for k in range(self.__n_x * self.__n_y):
                j = math.floor(k / self.__n_x)
                i = k - j * self.__n_x
                if i == 0:
                    self.__s[k, k] = 1
                    self.__s[k, k + 1] = -1
                elif i == self.__n_x - 1:
                    self.__s[k, k] = 1
                    self.__s[k, k - 1] = -1
                elif j == 0:
                    self.__s[k, k] = 1
                    self.__s[k, k + self.__n_x] = -1
                elif j == self.__n_y - 1:
                    self.__s[k, k] = 1
                    self.__s[k, k - self.__n_x] = -1
                else:
                    self.__s[k, k] = self.__parameter_c
                    self.__s[k, k + 1] = self.__s[k, k - 1] = -self.__parameter_d
                    self.__s[k, k + self.__n_x] = self.__s[k, k - self.__n_x] = -self.__parameter_e
        else:
            for k in range(self.__n_x * self.__n_y):
                j = math.floor(k / self.__n_x)
                i = k - j * self.__n_x
                if i == 0 or i == self.__n_x - 1 or j == 0 or j == self.__n_y - 1:
                    self.__s[k, k] = 1
                else:
                    self.__s[k, k] = self.__parameter_c
                    self.__s[k, k + 1] = self.__s[k, k - 1] = -self.__parameter_d
                    self.__s[k, k + self.__n_x] = self.__s[k, k - self.__n_x] = -self.__parameter_e

    def start_simulation(self):
        self.__lu, self.__piv = lu_factor(self.__s)

    def calculate_next_step(self, counter):
        for k in range(self.__n_x * self.__n_y):
            j = math.floor(k / self.__n_x)
            i = k - j * self.__n_x
            if i == 0 or i == self.__n_x - 1 or j == 0 or j == self.__n_y - 1:
                self.__b[k] = 0
            else:
                self.__b[k] = self.__parameter_f * self.__u_actual[k] - self.__parameter_g * self.__u_previous[k] + self.__parameter_h * (self.__u_previous[k + 1] + self.__u_previous[k - 1]) + self.__parameter_j * (self.__u_previous[k + self.__n_x] + self.__u_previous[k - self.__n_x]) + self.__parameter_k * (self.__u_actual[k + 1] + self.__u_actual[k - 1]) + self.__parameter_l * (self.__u_actual[k + self.__n_x] + self.__u_actual[k - self.__n_x])
                if self.__forcing_term_check:
                    self.__b[k] += self.__calculate_m(counter, k, self.__n_x, self.__delta_x, self.__delta_y, self.__delta_t, self.__parameter_a, self.__parameter_b, self.__beta, self.__x_f, self.__y_f, self.__t_f, self.__sigma)
        self.__u_next = lu_solve((self.__lu, self.__piv), self.__b)
        self.__u_previous = self.__u_actual
        self.__u_actual = self.__u_next

    def get_actual_step(self):
        return self.__u_actual

    def get_first_step(self, n_x, n_y, zero_diaphragm_displacement):
        delta_x = 1.0 / (n_x - 1)
        delta_y = 1.0 / (n_y - 1)
        u_first = np.zeros(n_x * n_y, dtype=np.float64)
        if not zero_diaphragm_displacement:
            for k in range(n_x * n_y):
                j = math.floor(k / n_x)
                i = k - j * n_x
                if not (i == 0 or i == n_x - 1 or j == 0 or j == n_y - 1):
                    u_first[k] = (i * delta_x - (i * delta_x) ** 2) * (j * delta_y - (j * delta_y) ** 2)

        return u_first

    def __delta_kronecker(self, x1, x2, delta):
        if math.fabs(x1 - x2) < delta * 0.4:
            return 1.0
        return 0.0

    def __forcing_term(self, x, x_f, y, y_f, t, t_f, delta_x, delta_y, sigma):
        return 100 * math.exp(- ((t - t_f) ** 2)/(2 * sigma ** 2)) * self.__delta_kronecker(x, x_f, delta_x) * self.__delta_kronecker(y, y_f, delta_y)

    def __calculate_m(self, n, k, __n_x, delta_x, delta_y, delta_t, parameter_a, parameter_b, beta, x_f, y_f, t_f, sigma):
        j = math.floor(k / __n_x)
        i = k - j * __n_x
        previous_extortion = self.__forcing_term(delta_x * i, x_f, delta_y * j, y_f, (n - 2) * delta_t, t_f, delta_x, delta_y, sigma)
        actual_extortion = self.__forcing_term(delta_x * i, x_f, delta_y * j, y_f, (n - 1) * delta_t, t_f, delta_x, delta_y, sigma)
        next_extortion = self.__forcing_term(delta_x * i, x_f, delta_y * j, y_f, n * delta_t, t_f, delta_x, delta_y, sigma)
        return delta_t ** 2 * (beta * next_extortion + parameter_a * previous_extortion + parameter_b * actual_extortion)
