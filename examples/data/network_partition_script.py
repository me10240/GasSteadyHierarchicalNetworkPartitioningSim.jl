import json
import logging
import networkx as nx
from networkx.algorithms.components import is_biconnected
log = logging.getLogger(__name__)
log.debug(nx.__version__)


def load_data_and_create_graph(filename):
    with open(filename, "r") as read_file:
        network_data = json.load(read_file)

    slack_nodes = []

    G = nx.Graph()
    G.clear()
    comp_list = ['pipes', 'short_pipes', 'compressors', 'valves', 'control_valves', 'resistors','loss_resistors']
    for comp in comp_list:
        if comp in network_data:
            for id in network_data[comp]:
                G.add_edge(network_data[comp][id]['fr_node'], network_data[comp][id]['to_node'])

    for node_id in network_data['nodes']:
        G.nodes[int(node_id)]["pos"] = (network_data['nodes'][node_id]['x_coord'], network_data['nodes'][node_id]['y_coord'])
        if network_data['nodes'][node_id]['slack_bool'] == 1:
            slack_nodes.append(int(node_id))
    slack_id = slack_nodes[0] #38, if multiple slacks, pick one of them 

    log.info("Network has {} nodes, {} edges, and the slack node(s) is/are {}".format(G.number_of_nodes(), G.number_of_edges(), slack_nodes))
    return G, slack_id




def create_list_of_articulation_points(G):

    if is_biconnected == True:
        log.warning("Network has no  articulation points")
        return []
    
    P = list(nx.articulation_points(G))
    def my_degree(n):
        return G.degree[n]
    P.sort(reverse=True, key=my_degree)
    log.info("Articulation points:{}".format(P))
    return P

def partition_graph_into_subgraphs_at_chosen_articulation_points(G, interface_nodes):
   
    if type(interface_nodes) is int:
        interface_node_list = [interface_nodes]
    else:
        interface_node_list = interface_nodes

    log.debug("interface_node_list:{}".format(interface_node_list))
    nbr_sets = []
    for i in interface_node_list:
        nbr_sets.append( list(G.neighbors(i)) )
    log.debug("nbr_sets:{}".format(nbr_sets))


    for i in range(len(interface_node_list)):
        nbrs = nbr_sets[i]
        node = interface_node_list[i]
        log.debug("node:{}, nbrs:{}".format(node, nbrs))
        for i in nbrs:
            if (node, i) in G.edges():
                G.remove_edge(node, i)

    S = [G.subgraph(c).copy() for c in sorted(nx.connected_components(G), key=len, reverse=False)]

    # assert exactly  len(interface_node_list) are singletons
    length_list = [bool(S[i].number_of_nodes() - 1)  for i in range(len(interface_node_list)+1) ]
    if any(length_list[:-1]) or not any(length_list):
        log.error("Partitioning issue! Check subnetwork sizes")
        exit()

    # remove first len(interface_node_list) subnetworks since they will be singleton interface nodes
    S = S[len(interface_node_list):]
        
    for ni, SG in enumerate(S):
        for i, intf_node in enumerate(interface_node_list): 
            nbrs = nbr_sets[i]
            for node in nbrs:
                if node in SG.nodes():
                    SG.add_edge(intf_node, node) 
    return S

def partition_graph(G):
    log.info("Nodes:{}, edges:{}".format(G.number_of_nodes(), G.number_of_edges()))
    P = create_list_of_articulation_points(G)
    if P == []:
        log.warning("Biconnected graph ! Could not partition.")
        return [], []
    interface_node_list = P[0:1]
    S  = partition_graph_into_subgraphs_at_chosen_articulation_points(G, interface_node_list)
    return S, interface_node_list

def put_slack_network_first(S, slack_id):
     # now knowing slack node, put subnetwork with slack node as first one
    for ni, SG in enumerate(S):
        if slack_id not in SG.nodes():
            log.debug("Slack node {} NOT in subnetwork {}+1 ".format(slack_id, ni))
            continue
        else:
            log.debug("Found slack node {} in subnetwork {}+1 ".format(slack_id, ni))
            if ni != 0:    
                S[0], S[ni] = S[ni], S[0] # swap 
            break
    return S

def write_partition_json_file_for_julia(filename, S, interface_node_list):
    partition_dict = {}
    # partition_dict["interface_node_list"] = interface_node_list
    for ni, SG in enumerate(S):
        partition_dict[ni+1] = list(SG.nodes()) #starting from 1 instead of 0
    log.info("Writing to json...")
    with open(filename, "w") as outfile:
        json.dump(partition_dict, outfile)
    log.info("Completed writing to json...")
    return partition_dict

def construct_block_cut_tree(p_dict, interface_node_list):

    G = nx.Graph()

    # G.add_nodes_from(interface_node_list)
    # G.add_nodes_from(list(range(1, len(p_dict)+1)))
    for i in range(1, len(p_dict)+1):
        for j in interface_node_list:
            if j in set(p_dict[i]):
                G.add_edge(i, j)
    
    
    import matplotlib.pyplot as plt
    plt.figure(figsize=(4, 4), dpi=200)
    color_list = ['seagreen' if node_name in interface_node_list else 'tan' for node_name in list(G.nodes)]
    size_list = [50 if node_name in interface_node_list else 200 for node_name in list(G.nodes)]
    font_size_list = [3 if node_name in interface_node_list else 6 for node_name in list(G.nodes)]


    pos = nx.spring_layout(G)
    nx.draw_networkx(G, pos=pos, with_labels=True, node_size= size_list, font_size=5, font_weight='bold', node_color=color_list, edge_color="black")
    plt.show()

    return


   



# show each partition with diff color at each stage ? possible ? interface nodes will be repeated in colouring
# import matplotlib.pyplot as plt
# m = len(S)
# plt.figure(figsize=(12, 12))
# plt.subplots_adjust(hspace=0.5)
# plt.suptitle("Partition", fontsize=18, y=0.95)
# for ni, SG in enumerate(S):
#     ax = plt.subplot(1, m, ni+1)
#     color_list = ['salmon' if node_name in interface_node_list else 'seagreen' if node_name == slack_id else 'bisque' for node_name in list(SG.nodes)]
#     nx.draw_spring(SG, with_labels=True, node_size= 600, font_size=16, font_weight='bold', node_color=color_list)

# plt.show()




def main():
    logging.basicConfig(format='%(asctime)s %(levelname)s--: %(message)s',
                        level=logging.INFO)
    
    filename = "GasLib-40-split/network.json"
    G, slack_id = load_data_and_create_graph(filename)
    S = [G]
    interface_node_list = []
    num_max = 20

    while True:
        remove_from_S = []
        add_to_S = []
        for index, SG in enumerate(S):
            if SG.number_of_nodes() > num_max:
                S_temp, interface_temp = partition_graph(SG)
                if S_temp == []:
                    continue
                remove_from_S.append(SG)
                add_to_S.extend(S_temp)
                interface_node_list.extend(interface_temp)
        if remove_from_S == []:
            break
        log.info("Partitioning ...")
        for item in remove_from_S:
            S.remove(item)
        S.extend(add_to_S)

    S = put_slack_network_first(S, slack_id)
    partition_data_file = "GasLib-40-split/partition-test-script.json"
    partition_dict = write_partition_json_file_for_julia(partition_data_file, S, interface_node_list)
    
    construct_block_cut_tree(partition_dict, interface_node_list)
    
    # import matplotlib.pyplot as plt
    # import networkx as nx
    # plt.figure(figsize=(6, 6), dpi=200)
    # color_list = ['red','orange','yellow','green','blue','purple', 'olive', 'cyan', 'pink']

    # pos = nx.get_node_attributes(G, "pos")
    # # log.debug("pos:{}".format(pos))
    # pos = nx.spring_layout(G, pos=pos, k=0.25, seed=10)
    # # color_list = ['seagreen' if node_name == slack_id else 'tan' for node_name in list(G.nodes)]
    # # nx.draw_networkx(SG, pos=pos, with_labels=True, node_size= 200, font_size=6, font_weight='bold', node_color=color_list[i], edge_color="tan")
    # # plt.show()      

    # for i, SG in enumerate(S):
    #     pos_dict = {}
    #     for j in SG.nodes():
    #         pos_dict[j]= pos[j]
    # # color_list = ['seagreen' if node_name == slack_id else 'tan' for node_name in list(G.nodes)]
    #     nx.draw_networkx(SG, pos=pos_dict, with_labels=True, node_size= 200, font_size=6, font_weight='bold', node_color=color_list[i], edge_color="tan")
    # plt.show()      

if __name__ == "__main__":
    main()