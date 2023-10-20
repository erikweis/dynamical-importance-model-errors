import json
from pathlib import Path

DIR = Path()
if not (DIR/'src').exists():
    DIR = Path('../..')
DATA_DIR = DIR / 'data'
EMPIRICAL_NETWORKS_DIR = DATA_DIR / 'empirical_networks'
EDGELIST_DIR = DATA_DIR / 'edgelists'

def create_edgelists():

    for path in EMPIRICAL_NETWORKS_DIR.iterdir():
        if path.name.startswith('.'):
            continue
        adj_list = json.load(open(path/'adjlist.txt'))
        edgelist_path = EDGELIST_DIR / rf'{path.name}.txt'
        with open(edgelist_path,'w') as f:
            for i,js in enumerate(adj_list):
                for j in js:
                    print(i,j,file=f)


if __name__ == "__main__":
    create_edgelists()

