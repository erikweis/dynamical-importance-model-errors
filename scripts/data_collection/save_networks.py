import graph_tool.all as gt
from pathlib import Path
import json

DIR = Path()
if not (DIR/'figures').exists():
    DIR = Path('../..')


def save_netzchleuder_graph(name,remove_self_loops = False):

    save_name = name.replace("/","_")
    print(name)
    g = gt.collection.ns[name]
    print("has self loops: ",sum(gt.label_self_loops(g)))
    if remove_self_loops:
        gt.remove_self_loops(g)
    
    network_dir = DIR / 'data/networks' / save_name
    network_dir.mkdir(exist_ok = True)

    # save adjlist
    adj_list = []
    for v in g.get_vertices():
        neighbors = [int(n) for n in g.get_all_neighbors(v)]
        adj_list.append(neighbors)
    with open(network_dir / 'adjlist.txt','w') as f:
        json.dump(adj_list, f)

    # save edgelist
    edge_list = [(int(edge[0]), int(edge[1])) for edge in g.get_edges()]
    with open(network_dir / 'edgelist.txt', 'w') as file:
        for edge in edge_list:
            print(*edge, file=file)


if __name__ == "__main__":

    networks = [
        "facebook_friends",
        "product_space/SITC",
        "faa_routes",
        "celegans_2019/male_gap_junction",
        "dolphins",
        "train_terrorists",
        "urban_streets/vienna",
        "student_cooperation",
        "polbooks",
        "malaria_genes/HVR_1",
        "lesmis",
        "karate/78",
        "jazz_collab",
        "euroroad",
        "crime"
    ]
    for network in networks:
        save_netzchleuder_graph(network, remove_self_loops=True)

    #g = gt.collection.ns['faa_routes']
    #print(gt.edge_reciprocity(g))
    #save_netzchleuder_graph("route_views_19990111",gt_name = "route_views/19990111", remove_self_loops = True)
    #save_netzchleuder_graph("physics_collab_pierreAuger",gt_name = "physics_collab/pierreAuger", remove_self_loops = False)
    #save_netzchleuder_graph("student_cooperation", remove_self_loops = False)
    #save_netzchleuder_graph("jazz_collab", remove_self_loops = False)
    #save_netzchleuder_graph("football_tsevans", remove_self_loops = False)
    #save_netzchleuder_graph("faa_routes")
    #save_netzchleuder_graph("facebook_friends")