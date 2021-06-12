"""
计算各种矩阵
"""
from math import sin, cos, pi
import numpy as np

# 构建节点功率方程
# 首先在Y矩阵中获取与其关联的节点，然后开始构建


#  PowerEquation类需要传入六个参数，已知功率，电压（模值）列向量，角度列向量，导纳矩阵，总节点数，PQ节点数m
class PowerEquation:
    def __init__(self, known_pq, voltage_diagonal, voltage_angle, admittance_matrix, shape, num_of_pq, DG, ub):
        self.known_pq = known_pq  # 已知的功率 [已知有功] + [已知无功]
        self.u = voltage_diagonal  # 节点电压幅值向量
        self.angle = voltage_angle  # 节点电压相角向量
        self.admittance = admittance_matrix  # 导纳矩阵
        self.shape = shape  # 总节点数量
        self.m = num_of_pq  # pq节点数量
        self.DG = DG
        self.ub = ub

    def P(self, *args):
        if not args:
            calculation_p = []
            for i in range(self.shape-1):  # 数组索引从零开始
                mutual_bus = [j for j, x in enumerate(self.admittance[i]) if x != 0]
                sum_p = 0
                for j in mutual_bus:
                    sum_p += self.u[j] * (self.admittance[i][j].real * cos(self.angle[i]-self.angle[j])
                                          + self.admittance[i][j].imag * sin(self.angle[i]-self.angle[j]))  # 有功功率计算方程
                sum_p *= self.u[i]
                calculation_p.append(sum_p)
            return calculation_p
        else:
            i = args[0]
            mutual_bus = [j for j, x in enumerate(self.admittance[i]) if x != 0]
            sum_p = 0
            for j in mutual_bus:
                sum_p += self.u[j] * (self.admittance[i][j].real * cos(self.angle[i] - self.angle[j])
                                      + self.admittance[i][j].imag * sin(self.angle[i] - self.angle[j]))  # 有功功率计算方程
            sum_p *= self.u[i]
            return sum_p

    def Q(self, *args):
        if not args:
            calculation_q = []
            for i in range(self.m):  # 数组索引从零开始
                mutual_bus = [j for j, x in enumerate(self.admittance[i]) if x != 0]
                sum_q = 0
                for j in mutual_bus:
                    sum_q += self.u[j] * (self.admittance[i][j].real * sin(self.angle[i] - self.angle[j])
                                          - self.admittance[i][j].imag * cos(self.angle[i] - self.angle[j]))  # 有功功率计算方程
                sum_q *= self.u[i]
                calculation_q.append(sum_q)
            return calculation_q
        else:
            i = args[0]
            mutual_bus = [j for j, x in enumerate(self.admittance[i]) if x != 0]
            sum_q = 0
            for j in mutual_bus:
                sum_q += self.u[j] * (self.admittance[i][j].real * sin(self.angle[i] - self.angle[j])
                                      - self.admittance[i][j].imag * cos(self.angle[i] - self.angle[j]))  # 有功功率计算方程
            sum_q *= self.u[i]
            return sum_q

    def get_delta_pq(self):
        delta_p = self.P()
        delta_q = self.Q()
        return self.known_pq - np.concatenate([delta_p, delta_q], axis=0)

    def convert_node(self, count):
        if self.DG:
            all_done = False
            for each in self.DG:
                if each[1] == 'PI':
                    print('{} \033[1;33mDG#{}-PI-convert to PQ\033[0m'.format(count, each[0]))
                    p = self.known_pq[each[0]] * 1e8
                    u = self.u[each[0]] * self.ub
                    i = each[2]
                    self.known_pq[self.shape - 1 + each[0]] = ((i**2 * u**2 - p ** 2) ** 0.5) / 1e8
                elif each[1] == 'PQV':
                    p = self.known_pq[each[0]] * 1e5
                    n = each[2]
                    u = self.u[each[0]] * 0.69 * 1000
                    q_export = n * 1000 * (u / 690)**2
                    # wind turbine connect to network through a transformer
                    # transformer params: 0.69kV/12.66kV, means Un = 690V
                    r = 5
                    x_t = 5
                    x_m = 2
                    s = r / (2 * p * (x_t ** 2)) * (u ** 2 - (u ** 4 - 4 * x_t ** 2 * p ** 2) ** 0.5)
                    q_import = p * (r ** 2 + x_t * (x_m + x_t) * s**2) / (r * x_m * s)

                    cos_phi_1 = p / (p**2 + (q_export - q_import)**2)**0.5
                    cos_phi_n = 0.95
                    if abs(cos_phi_1) < 0.95:
                        q_unit = 1000
                        q_export_change = p * (((1 / cos_phi_1) ** 2 - 1) ** 0.5 - ((1 / cos_phi_n) ** 2 - 1) ** 0.5)
                        if q_import > q_export:  # Need more PFC
                            add_n = int(q_export_change / q_unit) + 1
                            n += add_n
                            each[2] = n
                            q_export = n * q_unit * (u / 690)**2
                            self.known_pq[self.shape - 1 + each[0]] = (q_export - q_import)/1e8
                            # print('cosphi:{}'.format(p / (p**2 + (n * q_unit * (u / 690)**2 - q_import)**2)**0.5))
                            print('{} \033[1;33mDG#{}-PQV-{} sets of PFC were put in.'
                                  'Total sets:{}\033[0m'.format(count, each[0], add_n, n))
                        else:
                            sub_n = int(q_export_change / q_unit) + 1
                            n -= sub_n
                            each[2] = n
                            q_export = n * q_unit * (u / 690) ** 2
                            self.known_pq[self.shape - 1 + each[0]] = (q_export - q_import)/1e8
                            print('{} \033[1;33mDG#{}-PQV-{} sets of PFC were cut off.'
                                  'Total sets:{}\033[0m'.format(count, each[0], sub_n, n))
                    else:
                        print('{} \033[1;33mDG#{}-PQV-power factor:{}. No PFC changed.'
                              'Total sets:{}\033[0m'.format(count, each[0], round(cos_phi_1, 5), n))
                else:
                    if not all_done:
                        print('{} \033[1;36mAll DGs have converted successfully!\033[0m'.format(count))
                        all_done = True
        else:
            print('{} \033[1;34mNo DG Integrated\033[0m'.format(count))

    def jacobian_matrix_h(self):
        H = np.zeros((self.shape - 1, self.shape - 1), dtype=float)  # 确定雅可比矩阵的维度
        for i in range(self.shape-1):
            mutual_bus = [j for j, x in enumerate(self.admittance[i]) if x != 0 and j < self.shape-1]
            for j in mutual_bus:
                if i == j:
                    H[i][i] = (self.u[i] ** 2) * self.admittance[i][i].imag + self.Q(i)
                else:
                    H[i][j] = - self.u[i] * self.u[j] * (self.admittance[i][j].real * sin(self.angle[i] - self.angle[j])
                                      - self.admittance[i][j].imag * cos(self.angle[i] - self.angle[j]))
        return H

    def jacobian_matrix_n(self):
        N = np.zeros((self.shape-1, self.m), dtype=float)  # 确定雅可比矩阵的维度
        for i in range(self.shape-1):
            mutual_bus = [j for j, x in enumerate(self.admittance[i]) if x != 0 and j < self.m]
            for j in mutual_bus:
                if i == j:
                    N[i][i] = - (self.u[i] ** 2) * self.admittance[i][i].real - self.P(i)
                else:
                    N[i][j] = - self.u[i] * self.u[j] * (self.admittance[i][j].real * cos(self.angle[i]-self.angle[j])
                                      + self.admittance[i][j].imag * sin(self.angle[i]-self.angle[j]))
        return N

    def jacobian_matrix_m(self):
        M = np.zeros((self.m, self.shape-1), dtype=float)  # 确定雅可比矩阵的维度
        for i in range(self.m):
            mutual_bus = [j for j, x in enumerate(self.admittance[i]) if x != 0 and j < self.shape-1]
            for j in mutual_bus:
                if i == j:
                    M[i][i] = (self.u[i] ** 2) * self.admittance[i][i].real - self.P(i)
                else:
                    M[i][j] = self.u[i] * self.u[j] * (self.admittance[i][j].real * cos(self.angle[i]-self.angle[j])
                                      + self.admittance[i][j].imag * sin(self.angle[i]-self.angle[j]))
        return M

    def jacobian_matrix_l(self):
        L = np.zeros((self.m, self.m), dtype=float)  # 确定雅可比矩阵的维度
        for i in range(self.m):
            mutual_bus = [j for j, x in enumerate(self.admittance[i]) if x != 0 and j < self.m]
            for j in mutual_bus:
                if i == j:
                    L[i][i] = (self.u[i] ** 2) * self.admittance[i][i].imag - self.Q(i)
                else:
                    L[i][j] = - self.u[i] * self.u[j] * (self.admittance[i][j].real * sin(self.angle[i] - self.angle[j])
                                      - self.admittance[i][j].imag * cos(self.angle[i] - self.angle[j]))
        return L

    def get_jacobian_matrix(self):
        H = self.jacobian_matrix_h()
        N = self.jacobian_matrix_n()
        M = self.jacobian_matrix_m()
        L = self.jacobian_matrix_l()
        jacobian_matrix = - np.concatenate(
                    [np.concatenate([H, N], axis=1), np.concatenate([M, L], axis=1)], axis=0)
        return jacobian_matrix
