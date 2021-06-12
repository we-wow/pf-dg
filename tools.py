"""
Data为自定义的数据处理类，用于对输入数据进行处理，并获取节点导纳矩阵。
初始化时需要指定数据文件路径和节点导纳矩阵大小即节点数量即节点导纳矩阵大小即节点数量即:
Data(path=path, shape=shape)

PowerEquation为自定义潮流方程计算类
属性：
    各个电压节点的幅值：u
    角度：angle
    导纳矩阵：admittance
方法:
    计算△PQ的方法：get_delta_pq()
    获取雅克比矩阵的方法:get_jacobian_matrix()

初始值设置为平启动/平直电压法
"""
import math
from read import Data
from matrix import PowerEquation
from functools import wraps
import numpy as np
import time
import copy


def fn_timer(function):
    @wraps(function)
    def function_timer(*args, **kwargs):
        t0 = time.time()
        result = function(*args, **kwargs)
        t1 = time.time()
        print("Total time running %s: %s seconds" %
              (function.__name__, str(t1 - t0))
              )
        return result
    return function_timer


def data_initial(branch_path, power_path, u_base=12660):
    data = Data(branch_file_path=branch_path, node_file_path=power_path, u_base=u_base)
    return data


def network_initial(data):
    net_work = PowerEquation(
        known_pq=data.known_pq,
        voltage_diagonal=data.voltage,
        voltage_angle=data.angle,
        admittance_matrix=data.admittance_matrix,  # 导纳矩阵
        shape=data.shape,  # 节点数，导纳矩阵维数
        num_of_pq=data.num_of_pq,  # PQ节点数
        DG=data.DG,
        ub=data.u_base
    )
    return net_work


def data_return(shape, no2no, net_work):
    voltage = {}  # 返回节点电压
    angle = {}  # 返回节点相角
    active_power = {}  # 返回节点有功
    reactive_power = {}  # 返回节点无功
    get_p = []
    get_q = []
    for i in range(shape):
        get_p.append(net_work.P(i))
        get_q.append(net_work.Q(i))
    for i in range(0, shape):
        index = no2no[i]  # 获取节点对应运算编号
        radian2angle = []
        for radian in net_work.angle:
            radian2angle.append(radian / 2 / math.pi * 360)  # 弧度转换为角度
        voltage[i] = net_work.u[index]
        active_power[i] = get_p[index]
        reactive_power[i] = get_q[index]
        angle[i] = radian2angle[index]
    return {'voltage': voltage,
            'angle': angle,
            'active_power': active_power,
            'reactive_power': reactive_power,
            }


@fn_timer
def newton_method(data, net_work):
    # 计算，初始化网络方程类
    x = data.angle[:data.shape-1] + data.voltage[:data.num_of_pq]
    count = 0
    while count < 50:
        # X(k) <---- X(k-1)
        # 获取△PQ和雅克比矩阵
        delta_pq = net_work.get_delta_pq()
        if max(abs(delta_pq)) <= 1e-9:
            break  # 设置退出条件
        jacobian_matrix = net_work.get_jacobian_matrix()
        # 计算△X
        delta_x = np.linalg.solve(jacobian_matrix, delta_pq)
        delta_x[data.shape-1:] *= net_work.u[:data.num_of_pq]
        # 计算X(k+1)
        x += delta_x
        # 更新数据，角度和电压幅值
        net_work.angle[:data.shape-1] = x[:data.shape-1]
        net_work.u[:data.num_of_pq] = x[data.shape-1:]
        count += 1
        net_work.convert_node(count)
    if count == 50:
        print("Data divergence...")
    data4return = data_return(shape=data.shape, no2no=data.o2n, net_work=net_work)
    return data4return, count, data


