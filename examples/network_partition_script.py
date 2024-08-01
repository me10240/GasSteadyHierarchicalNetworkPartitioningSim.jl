import json
import logging
import networkx as nx
from networkx.algorithms.components import is_biconnected
log = logging.getLogger(__name__)
log.debug(nx.__version__)


def load_data_and_create_graph(filename, plotting_flag=False):
    with open(filename, "r") as read_file:
        network_data = json.load(read_file)

    slack_nodes = []

    G = nx.Graph()
    G.clear()
    comp_list = ['pipes', 'short_pipes', 'compressors', 'valves', 'control_valves', 'resistors','loss_resistors']
    for comp in comp_list:
        if comp in network_data:
            for id in network_data[comp]:
                if "fr_node" in network_data[comp][id].keys():
                    G.add_edge(network_data[comp][id]['fr_node'], network_data[comp][id]['to_node'])
                else:
                    G.add_edge(network_data[comp][id]['from_node'], network_data[comp][id]['to_node'])

    for node_id in network_data['nodes']:
        G.nodes[int(node_id)]["pos"] = (network_data['nodes'][node_id]['x_coord'], network_data['nodes'][node_id]['y_coord'])
        if network_data['nodes'][node_id]['slack_bool'] == 1:
            slack_nodes.append(int(node_id))

    log.info("Network has {} nodes, {} edges, and the slack node(s) is/are {}".format(G.number_of_nodes(), G.number_of_edges(), slack_nodes))

    if plotting_flag:
        pos = nx.get_node_attributes(G, "pos")
        import matplotlib.pyplot as plt
        plt.figure(figsize=(4.5, 4.5), dpi=200)
        color_list = ["seagreen" if node_name in slack_nodes  else 'darkorange' for node_name in list(G.nodes)]
        nx.draw_networkx(G, pos=pos, with_labels=False, node_size= 5, font_size=5, font_weight='bold', node_color=color_list, edge_color="black")
        plt.show()
    return G, slack_nodes


def create_list_of_articulation_points(G):

    P = list(nx.articulation_points(G))
    def my_degree(n):
        return G.degree[n]
    P.sort(reverse=True, key=my_degree)
    log.debug("Articulation points:{}".format(P))
    return P

def partition_graph_into_subgraphs_at_chosen_articulation_points(G, interface_nodes, allow_singletons_flag=True):
   
    if type(interface_nodes) == int:
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

    removal_list = []
    for SG in S:
        for node in interface_node_list:
            if node in SG.nodes():
                removal_list.append(SG)
                break
    for SG in removal_list:
        S.remove(SG)

    if not allow_singletons_flag:
        # check if any singletons remain
        if S[0].number_of_nodes() == 1:
            log.error("Partitioning issue! Singletons found !")
            exit()

        
    for ni, SG in enumerate(S):
        for i, intf_node in enumerate(interface_node_list): 
            nbrs = nbr_sets[i]
            for node in nbrs:
                if node in SG.nodes():
                    SG.add_edge(intf_node, node) 
    return S

def partition_graph(G, slack_nodes, allow_slack_node_partitioning=True):

    log.debug("Nodes:{}, edges:{}".format(G.number_of_nodes(), G.number_of_edges()))
    P = create_list_of_articulation_points(G)

    if P == []:
        log.warning("Found biconnected network with {} nodes !".format(G.number_of_nodes()))
        return [], []

    if not allow_slack_node_partitioning:
        P = list(set(P) - set(slack_nodes))
        if P == []:
            log.warning("Non-slack articulation points not found! Could not partition.")
            return [], []

    l = len(set(G.nodes()).intersection(set(slack_nodes)))
    for i in range(len(P)):
        interface_node = P[i]
        S  = partition_graph_into_subgraphs_at_chosen_articulation_points(G, interface_node)
        var =[len( set(slack_nodes).intersection(set(SG.nodes())) ) for SG in S]
        if max(var) == l:
            return S, [interface_node]
    
    return [], []
    

def put_slack_networks_first(S, slack_nodes):

    def slack_ordering(SG, slack_nodes):
        return len(set(SG.nodes()).intersection(set(slack_nodes)))
        
    S.sort(reverse=True, key=lambda SG: slack_ordering(SG, slack_nodes))

    slack_network_ids = []
    # identify slack networks
    for slack_id in slack_nodes:
        for ni, SG in enumerate(S):
            if slack_id in SG.nodes() and ni not in slack_network_ids:
                slack_network_ids.append(ni)
                log.debug("Found slack node {} in subnetwork {}+1 ".format(slack_id, ni))
    
    if len(slack_network_ids) == 1:
        log.info("All slack nodes within single subnetwork")
    else:
        log.warning("There are {} slack subnetworks".format(len(slack_network_ids)))

    return S, slack_network_ids

def write_partition_json_file_for_julia(filename, S, interface_node_list, slack_network_ids, slack_nodes):
    partition_dict = {}
    partition_dict["num_partitions"] = len(S)
    partition_dict["interface_nodes"] = interface_node_list
    partition_dict["slack_network_ids"] = [i+1 for i in slack_network_ids]
    partition_dict["slack_nodes"] = slack_nodes
    for ni, SG in enumerate(S):
        partition_dict[ni+1] = list(SG.nodes()) #starting from 1 instead of 0
    log.info("Writing to json...")
    with open(filename, "w") as outfile:
        json.dump(partition_dict, outfile)
    log.info("Completed writing to json...")
    return partition_dict

def construct_block_cut_tree(p_dict, plotting_flag=False, full_network=None):

    G = nx.Graph()
    num_partitions = p_dict["num_partitions"]
    interface_node_list = p_dict["interface_nodes"]
    slack_network_ids = p_dict["slack_network_ids"]

    for i in range(1, num_partitions+1):
        for j in interface_node_list:
            if j in set(p_dict[i]):
                node = "N-{}".format(i)
                G.add_edge(node, j)
    tree_status = nx.is_tree(G)
    log.info("Block cut graph is a tree: {}".format(tree_status))

    if plotting_flag:
        if full_network:
            def compute_partition_centroid(G0, p_dict, i):
                centroid_x = 0.0
                centroid_y = 0.0
                num_nodes = len(p_dict[i])
                for node in p_dict[i]:
                    centroid_x += G0.nodes[node]["pos"][0]
                    centroid_y += G0.nodes[node]["pos"][1]
                return (centroid_x/num_nodes, centroid_y/num_nodes)
            
            for i in range(1, num_partitions+1):
                G.nodes[f"N-{i}"]["pos"] = compute_partition_centroid(full_network, p_dict, i)
            for i in interface_node_list:
                G.nodes[i]["pos"] = full_network.nodes[i]["pos"]
            pos = nx.get_node_attributes(G, "pos")
        else:
            pos = nx.spring_layout(G)

        slack_network_names = [f"N-{i}" for i in slack_network_ids]
        import matplotlib.pyplot as plt
        plt.figure(figsize=(4.5, 4.5), dpi=200)
        color_list = ['tan' if node_name in interface_node_list else "seagreen" if node_name in slack_network_names  else 'darkorange' for node_name in list(G.nodes)]
        size_list = [50 if node_name in interface_node_list else 150 for node_name in list(G.nodes)]
        font_size_list = [2 if node_name in interface_node_list else 4 for node_name in list(G.nodes)]

        nx.draw_networkx(G, pos=pos, with_labels=True, node_size= size_list, font_size=5, font_weight='bold', node_color=color_list, edge_color="black")
        plt.show()

    return

def partition_given_network(dirname, allow_slack_node_partitioning=True, num_max=20, round_max=3, plotting_flag=False):
    filename = dirname + "network.json"
    G, slack_nodes = load_data_and_create_graph(filename, plotting_flag=plotting_flag)
    S = [G]
    interface_node_list = []
    nonseparable_subnetworks = []

    if num_max >= G.number_of_nodes():
        log.error("max partition size must be greater than network size!")
        exit()

    partitioning_round = 0
    while True:
        remove_nonseparable_from_S = []
        for index, SG in enumerate(S):
            if nx.is_biconnected(SG) == True:
                remove_nonseparable_from_S.append(SG)
        for item in remove_nonseparable_from_S:
            S.remove(item)
            nonseparable_subnetworks.append(item)

        if S == []:
            log.info("Cannot partition further, no subnetworks left")
            break

        if partitioning_round == round_max:
            log.info("Completed {} rounds of partitioning".format(round_max))
            break

        log.debug("Partitioning... round {}".format(partitioning_round))
        remove_from_S = []
        add_to_S = []

        def my_size(SG):
            return SG.number_of_nodes()

        S.sort(reverse=True, key=my_size)

        if S[0].number_of_nodes() < num_max:
            log.info("No need to partition block. Partitioning concludes after {} rounds ".format(partitioning_round))
            break
        else:
            S_temp, interface_temp = partition_graph(S[0], slack_nodes, allow_slack_node_partitioning=allow_slack_node_partitioning)

        if S_temp == []:
            nonseparable_subnetworks.append(S[0])
        else:
            S.extend(S_temp)
            interface_node_list.extend(interface_temp)

        S.remove(S[0])
        partitioning_round += 1
    S.extend(nonseparable_subnetworks)
    S, slack_network_ids = put_slack_networks_first(S, slack_nodes)

    partition_sizes = [SG.number_of_nodes()   for SG in S ]
    if max(partition_sizes) <= num_max:
        log.info("Partitioning achieved desired size")
    else:
        log.info("Partition has nonseparable  block larger than desired size")

    log.info("Size of largest partition is {}".format(max(partition_sizes)))

    if len(S) >= 2:
        partition_data_file = dirname + "partition-test-script.json"
        partition_dict = write_partition_json_file_for_julia(partition_data_file, S, interface_node_list, slack_network_ids, slack_nodes)
        construct_block_cut_tree(partition_dict, plotting_flag=plotting_flag, full_network=None)
    else:
        log.info("Could not partition network")
    return

def main():
    logging.basicConfig(format='%(asctime)s %(levelname)s--: %(message)s',
                        level=logging.INFO)
    logging.getLogger('matplotlib.font_manager').disabled = True
    
    dirname = "Texas7k_Gas/"
    partition_given_network(dirname, allow_slack_node_partitioning = False, num_max=2, round_max=100, plotting_flag=True)

def run_script(dirname, allow_slack_node_partitioning = False, num_max=2, round_max=1, plotting_flag=True):
    logging.basicConfig(format='%(asctime)s %(levelname)s--: %(message)s',
                        level=logging.INFO)
    logging.getLogger('matplotlib.font_manager').disabled = True
    
    # dirname = "8-node/"
    partition_given_network(dirname, allow_slack_node_partitioning = allow_slack_node_partitioning, num_max=num_max, round_max=round_max, plotting_flag=plotting_flag)

    
if __name__ == "__main__":
    main()