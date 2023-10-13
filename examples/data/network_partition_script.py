import networkx as nx
print(nx.__version__)

import json
with open("GasLib-40-split/network.json", "r") as read_file:
    network_data = json.load(read_file)

slack_nodes = []

for node_id in network_data['nodes']:
    if network_data['nodes'][node_id]['slack_bool'] == 1:
        slack_nodes.append(int(node_id))
print("The slack node(s) is/are ", [int(s) for s in slack_nodes])
slack_id = slack_nodes[0] #38, if multiple slacks, pick one of them 


G = nx.Graph()
G.clear()
for pipe_id in network_data['pipes']:
    G.add_edge(network_data['pipes'][pipe_id]['fr_node'], network_data['pipes'][pipe_id]['to_node'])
for comp_id in network_data['compressors']:
    G.add_edge(network_data['compressors'][comp_id]['fr_node'], network_data['compressors'][comp_id]['to_node'])

print(G.number_of_nodes(), G.number_of_edges())


max_degree_index = 1
for i in range(1, G.number_of_nodes() + 1):
    if G.degree[i] > max_degree_index:
        max_degree_index = G.degree[i]
    else:
        continue
print(max_degree_index)

max_degree_nodes = []
for i in range(1, G.number_of_nodes() + 1):
    if G.degree[i] == max_degree_index:
        max_degree_nodes.append(i)
    else:
        continue
print(max_degree_nodes)


interface_node_list = [21, 27, 32]
nbr_sets = []
for i in interface_node_list:
    nbr_sets.append( list(G.neighbors(i)) )
print(nbr_sets)


for i in range(len(interface_node_list)):
    nbrs = nbr_sets[i]
    node = interface_node_list[i]
    print(node, nbrs)
    for i in nbrs:
        if (node, i) in G.edges():
            G.remove_edge(node, i)

len_subnetworks_array = [len(c) for c in sorted(nx.connected_components(G), key=len, reverse=True)]
print(len_subnetworks_array)

num_isolated_nodes = 0
for num in len_subnetworks_array:
    if num == 1:
        num_isolated_nodes += 1

if num_isolated_nodes != len(interface_node_list):
    print("Number of isolated node subnetworks is {}, expect/want  {}".format(num_isolated_nodes, len(interface_node_list)))



S = [G.subgraph(c).copy() for c in sorted(nx.connected_components(G), key=len, reverse=False)]

# remove first len(interface_node_list) subnetworks since they will be singleton interface nodes
S = S[len(interface_node_list):]
for SG in S:
    print(SG.nodes(), SG.edges())
    

# # now knowing slack node, put subnetwork with slack node as first one

for ni, SG in enumerate(S):
    if slack_id not in SG.nodes():
        print("Slack node {} NOT in subnetwork {}+1 ".format(slack_id, ni))
        continue
    else:
        print("Found slack node {} in subnetwork {}+1 ".format(slack_id, ni))
        
    if ni == 0:
        break 
        
    SG_temp = S[0]
    S[0] = SG
    S[ni] = SG_temp
    break
        
partition_dict = {}
for ni, SG in enumerate(S):
    for i in range(len(interface_node_list)):
        intf_node = interface_node_list[i]
        nbrs = nbr_sets[i]
        for node in nbrs:
            if node in SG.nodes():
                SG.add_edge(intf_node, node) 
    partition_dict[ni+1] = list(SG.nodes())

import matplotlib.pyplot as plt
m = len(S)
plt.figure(figsize=(12, 12))
plt.subplots_adjust(hspace=0.5)
plt.suptitle("Partition", fontsize=18, y=0.95)
for ni, SG in enumerate(S):
    ax = plt.subplot(1, m, ni+1)
    color_list = ['salmon' if node_name in interface_node_list else 'seagreen' if node_name == slack_id else 'bisque' for node_name in list(SG.nodes)]
    nx.draw_spring(SG, with_labels=True, node_size= 600, font_size=16, font_weight='bold', node_color=color_list)

plt.show()


with open("GasLib-40-split/partition-dummy.json", "w") as outfile:
    json.dump(partition_dict, outfile)