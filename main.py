from tools import *
data = data_initial(branch_path='./data/branch.csv', power_path='./data/node.csv', u_base=12660)

# [father=0, no=0, r=0.1, x=0.1, p=0.0001, q=0.0001, v=1.0, node_type='PQ', more=0.0]
# this is scenario 5's data.
data.add_node([19, 33, 0.1, 0.1, 0.01, 0.0, 1.0, 'PQV', 10])
data.add_node([17, 33, 0.1, 0.1, 0.004, 0.0, 1.0, 'PI', 50])
data.add_node([29, 33, 0.1, 0.1, 0.01, 0.0, 1.0, 'PV', 0])
data.add_node([24, 36, 0.1, 0.1, 0.005, 0.005, 1.0, 'DGPQ', 0])

network = network_initial(data)
result, times, source = newton_method(data, network)
loss = data.calculate_loss(result['voltage'], result['angle'])
print('Total Loss:{}MVA'.format(loss))
data.draw_network(result)
data.display_result(result)
data.write2csv(data=result, file_path='IEEE33_result_case_5.csv')
