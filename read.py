import csv
import numpy as np
import copy
import math

class Node:
    def __init__(self, no=0, voltage=1.0, father=0):
        self.no = no
        self.voltage = voltage  # random.uniform(0.8, 1)
        self.father = father
        self.sons = []
        self.check_sons = copy.deepcopy(self.sons)
        self.check_father = father


class Data:
    def __init__(self, branch_file_path, node_file_path, u_base=12660):
        self.u_base = u_base
        self.branch_file_path = branch_file_path
        self.node_file_path = node_file_path
        self.branch_read_from_file = self.read_branch_file()
        self.node_read_from_file = self.read_node_file()
        self.shape = len(self.node_read_from_file)
        self.n2o, self.o2n, self.num_of_pq = self.rename_bus()
        self.admittance_matrix = self.get_admittance_matrix()
        self.known_pq, self.voltage, self.angle, self.DG= self.get_power_data()

    def read_branch_file(self):
        with open(self.branch_file_path, 'r', encoding='utf-8-sig') as data:
            data = csv.reader(data)
            temp_data = []
            for each in data:
                temp_data.append([int(each[0]), int(each[1]), float(each[2]), float(each[3])])
        return temp_data

    def read_node_file(self):
        with open(self.node_file_path, 'r', encoding='utf-8-sig') as data:
            data = csv.reader(data)
            temp_data = []
            for each in data:
                temp_data.append([int(each[0]), float(each[1]), float(each[2]),
                                  each[3], float(each[4]), eval(each[5])])
        return temp_data

    def rename_bus(self):
        pq_bus = []
        pv_bus = []
        slack_bus = 0
        branch_data = self.branch_read_from_file
        node_data = self.node_read_from_file
        for each in node_data:
            i = each[0]
            if each[3] == 'SL':  # 参考节点
                slack_bus = i
            elif each[3] == 'PV':  # PV节点
                pv_bus.append(i)
            else:  # PQ节点
                pq_bus.append(i)
        n2o = {}
        o2n = {}
        search_list = pq_bus + pv_bus + [slack_bus]
        value = 0
        for i in search_list:
            o2n[i] = value
            n2o[value] = i
            value += 1
        for node in node_data:
            node[0] = o2n[node[0]]
        node_data.sort(key=lambda x: (x[0]))
        for branch in branch_data:
            branch[0] = o2n[branch[0]]
            branch[1] = o2n[branch[1]]
        return n2o, o2n, len(pq_bus),

    def get_admittance_matrix(self):
        data_array = np.zeros((self.shape, self.shape), dtype=complex)
        branch_data = self.branch_read_from_file
        for each in branch_data:
            i = each[0]
            j = each[1]
            r = each[2]
            x = each[3]
            Yij = 1 / complex(r, x)
            data_array[i][i] += Yij
            data_array[i][j] -= Yij
            data_array[j][i] -= Yij
            data_array[j][j] += Yij
        return data_array

    def get_power_data(self):
        node_data = self.node_read_from_file
        known_pq = [0.0] * (self.shape-1 + self.num_of_pq)
        voltage = [1.0] * self.shape
        angle = [0.0] * self.shape
        DG = []
        for node in node_data:
            if node[3] != 'SL':
                known_pq[node[0]] = node[1]  # read the active power
            if node[3] == 'PQ':
                known_pq[self.shape - 1 + node[0]] = node[2]  # read the reactive power
            elif node[3] == 'DGPQ':
                known_pq[self.shape - 1 + node[0]] = node[2]  # read the reactive power
                DG.append([node[0], 'DGPQ'])
            elif node[3] == 'PI':
                known_pq[self.shape - 1 + node[0]] = \
                    ((node[5]**2 * (node[4]*self.u_base)**2 - (node[1]*1e8)**2)**0.5) / 1e8  # read the reactive power
                DG.append([node[0], 'PI', node[5]])
            elif node[3] == 'PQV':
                known_pq[self.shape - 1 + node[0]] = node[2]  # read the reactive power
                DG.append([node[0], 'PQV', node[5]])
            elif node[3] == 'PV':
                voltage[node[0]] = node[4]  # read the voltage
                DG.append([node[0], 'PV'])
        return known_pq, voltage, angle, DG

    def add_node(self, *args):
        """
        eg:[father=0, no=0, r=0.1, x=0.1, p=0.0001, q=0.0001, v=1.0, node_type='PQ', more=0.0]
        [node1],[node2],[node3]......
        The type of each params should be [int int float float float float float float str float]
        For example, it should be V=1.0 rather than V=1.
        """
        if args:
            branch_data = self.branch_read_from_file
            node_data = self.node_read_from_file
            for branch in branch_data:
                branch[0] = self.n2o[branch[0]]
                branch[1] = self.n2o[branch[1]]
            for node in node_data:
                node[0] = self.n2o[node[0]]
            for each in args:
                if (len(each) == 9 and type(each[0]) == int and type(each[1]) == int and type(each[2]) == float
                        and type(each[3]) == float and type(each[4]) == float and type(each[5]) == float
                        and type(each[6]) == float and type(each[7]) == str):
                    self.branch_read_from_file.append([each[0], each[1], each[2], each[3]])
                    self.node_read_from_file.append([each[1], each[4], each[5], each[7], each[6], each[8]])
                    print('\033[1;33mDG#{} ({}) integrated to node#{}\033[0m'.format(each[1], each[7], each[0]))
                else:
                    raise ValueError('Unexpected Type')
            self.shape = len(self.node_read_from_file)
            self.n2o, self.o2n, self.num_of_pq = self.rename_bus()
            self.admittance_matrix = self.get_admittance_matrix()
            self.known_pq, self.voltage, self.angle, self.DG = self.get_power_data()

        else:
            Warning("Nothing Input")

    def draw_network(self, result={}):
        import networkx as nx
        import matplotlib.pyplot as plt
        import matplotlib
        voltage = [0.5]*(self.shape+1)
        angle = []
        p = []
        q = []
        voltage_given = False
        if result:
            voltage = result['voltage']
            angle = result['angle']
            p = result['active_power']
            q = result['reactive_power']
            voltage_given = True
        node = {}
        for line in self.node_read_from_file:
            n = Node(no=self.n2o[line[0]], voltage=voltage[self.n2o[line[0]]])
            node[self.n2o[line[0]]] = n
        for line in self.branch_read_from_file:
            node[self.n2o[line[1]]].father = self.n2o[line[0]]
            node[self.n2o[line[0]]].sons.append(self.n2o[line[1]])

        # to draw a picture of the structure of this network
        for each_node in node.values():
            each_node.check_father = copy.deepcopy(each_node.father)
            each_node.check_sons = copy.deepcopy(each_node.sons)
        grid = nx.DiGraph()
        index = list(range(1, len(node), 1))
        node[0].check_father = 'done'
        pos = {}

        axisUnit = 50

        x = axisUnit
        y = 0
        flag = 0
        pos[0] = (0, 0)
        while index:
            searchList = []
            fatherList = []
            for i in index:
                if node[node[i].father].check_father == 'done':
                    searchList.append(i)
                    fatherList.append(node[i].father)
            if len(node[fatherList[0]].sons) > 1:
                tempList = copy.deepcopy(searchList)
                searchList = copy.deepcopy(sorted(node[fatherList[0]].sons))
                tempList = sorted(list(set(searchList) ^ set(tempList)),
                                  reverse=False)  # Need sorted by his father's location
                sortList = []
                for each in tempList:
                    sortList.append((each, pos[node[each].father][1]))
                sortList.sort(key=lambda x: x[1])
                for n in range(len(sortList)):
                    tempList[n] = sortList[n][0]
                searchList = searchList + tempList

            locationOccupyed = []
            for i in searchList:
                if node[node[i].father].voltage >= node[i].voltage:
                    grid.add_edge(node[i].father, i)
                else:
                    grid.add_edge(i, node[i].father)
                if pos[node[i].father][1] > y:
                    y = pos[node[i].father][1]

                if y in locationOccupyed:
                    y += axisUnit
                    while (y in locationOccupyed):
                        y += axisUnit

                pos[i] = (x, y)
                locationOccupyed.append(y)
                node[i].check_father = 'done'
                index.remove(i)
                flag += 1
                y = flag * axisUnit
            for i in searchList:
                if pos[node[i].father][1] < pos[i][1]:
                    temp_y = copy.deepcopy(pos[node[i].father][1])
                    while (temp_y in locationOccupyed and temp_y != pos[i][1]):
                        temp_y += axisUnit
                    locationOccupyed.remove(pos[i][1])
                    locationOccupyed.append(temp_y)
                    pos[i] = (x, temp_y)
            x += axisUnit
            y = 0
            flag = 0
        pos[0] = (0, 0)
        cdict = {'red': ((0.0, 1.0, 1.0),
                         (1.0, 0.0, 0.0)),

                 'green': ((0.0, 0.0, 0.0),
                           (1.0, 1.0, 1.0)),

                 'blue': ((0.0, 0.0, 0.0),
                          (1.0, 0.0, 0.0))}
        blue_red2 = matplotlib.colors.LinearSegmentedColormap('BlueRed2', cdict)
        plt.register_cmap(name='blue_red', cmap=blue_red2)
        cmap = plt.get_cmap('blue_red')  # 'BuPu'
        plt.subplot(211)
        alpha = []
        for i in grid.nodes:
            alpha.append(node[i].voltage)
        nx.draw(grid, pos, node_color=alpha, node_size=400, with_labels=True, vmin=min(alpha), vmax=max(alpha),
                cmap=cmap)
        plt.axis('equal')
        real_voltage = [0] * len(node)
        for key, each in node.items():
            real_voltage[key] = each.voltage
        max_index = [i for i in range(len(real_voltage)) if abs(real_voltage[i] - max(real_voltage)) < 1e-6]
        min_index = [i for i in range(len(real_voltage)) if real_voltage[i] == min(real_voltage)]
        if not voltage_given:
            title = '{}nodes power system networks'.format(self.shape)
        else:
            title = '{}nodes power system networks voltage loss\n\n ' \
                'Green represents the max voltage of all. Red means it has the lowest voltage of any node\n ' \
                'For this system green to red is {}kV{} to {}kV{}'.format(
                len(node), round(max(real_voltage)*self.u_base/1000, 5), max_index,
                round(min(real_voltage)*self.u_base/1000, 5), min_index)
        plt.title(title)

        if voltage_given and self.DG:
            data = []
            for dg in self.DG:
                no = self.n2o[dg[0]]
                dg_type = dg[1]
                u = round(voltage[no] * 12.66, 5)
                a = round(angle[no], 5)
                dg_p = round(p[no] * 100, 5)
                dg_q = round(q[no] * 100, 5)
                data.append([no, dg_type, u, a, dg_p, dg_q])
            my_table = plt.table(
                cellText=data,
                colColours=['red', 'blue', 'red', 'blue', 'red', 'blue'],
                colLabels=('No', 'Type', 'Voltage/kV', 'Angle/°', 'P/MVA', 'Q/MVA'),
                loc='bottom',
            )
            my_table.auto_set_font_size(False)
            my_table.set_fontsize(10)
        plt.show()

    def calculate_loss(self, voltage, angle):
        loss = 0
        for each in self.read_branch_file():
            u1 = voltage[each[0]]
            a1 = angle[each[0]] / 180 * math.pi
            u2 = voltage[each[1]]
            a2 = angle[each[1]] / 180 * math.pi
            u1 = u1 * complex(math.cos(a1), math.sin(a1))
            u2 = u2 * complex(math.cos(a2), math.sin(a2))
            z = complex(each[2], each[3])
            branch_loss = (u1 - u2) * ((u1 - u2) / z).conjugate()
            # if branch_loss.real < 0 and branch_loss.imag < 0:
            #     branch_loss = -branch_loss
            loss += branch_loss
        return loss * 100

    def write2csv(self, data, file_path='./result.csv'):
        import csv
        u = data['voltage']
        a = data['angle']
        p = data['active_power']
        q = data['reactive_power']
        data2write = []
        for i in u.keys():
            u_i = round(u[i] * self.u_base/1000, 6)
            a_i = round(a[i], 6)
            p_i = round(p[i] * 100, 6)
            q_i = round(q[i] * 100, 6)
            data2write.append([i, u_i, a_i, p_i, q_i])
        header = ['No', 'V/kV', 'Angle/°', 'P/MVA', 'Q/MVA']
        with open(file_path, 'w', newline='') as f:
            f_csv = csv.writer(f)
            f_csv.writerow(header)
            f_csv.writerows(data2write)
        print('\033[1;34mData has been successfully written to {}\033[0m'.format(file_path))

    def display_result(self, data):
        import prettytable
        u = data['voltage']
        a = data['angle']
        p = data['active_power']
        q = data['reactive_power']
        data2write = []
        for i in u.keys():
            u_i = round(u[i] * self.u_base / 1000, 6)
            a_i = round(a[i], 6)
            p_i = round(p[i] * 100, 6)
            q_i = round(q[i] * 100, 6)
            data2write.append([i, u_i, a_i, p_i, q_i])
        header = ['No', 'V/kV', 'Angle/°', 'P/MVA', 'Q/MVA']
        tb = prettytable.PrettyTable()
        tb.field_names = header
        for each in data2write:
            tb.add_row(each)
        print(tb)



if __name__ == '__main__':
    data = Data(branch_file_path='./data/branch.csv', node_file_path='./data/node.csv')
    data.draw_network()