import matplotlib.pyplot as plt
import numpy as np
import random
import sys
import argparse

def read_tsp(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    coords = []
    reading_nodes = False
    for line in lines:
        if line.startswith("NODE_COORD_SECTION"):
            reading_nodes = True
            continue
        if reading_nodes:
            if line.startswith("EOF"):
                break
            parts = line.strip().split()
            if len(parts) >= 3:
                x = float(parts[1])
                y = float(parts[2])
                coords.append((x, y))
    return np.array(coords)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("tsp_file", help="Path to .tsp file")
    parser.add_argument("tsp_path", help="String of nodes of TSP path", default=None, nargs='?')
    args = parser.parse_args()

    cities = read_tsp(args.tsp_file)
    n_cities = len(cities)
    
    if args.tsp_path:
        # path is string "[1 2 3 ...]"
        path = [int(x)-1 for x in args.tsp_path.strip("[]").split(" ")]
    else:
        path = list(range(n_cities))
        random.shuffle(path)

    plt.figure(figsize=(8, 6))
    plt.scatter(cities[:, 0], cities[:, 1], color='blue', zorder=5)
    for i, (x, y) in enumerate(cities):
        plt.text(x, y, str(i+1), fontsize=10, ha='right', va='bottom')
    for i in range(len(path) - 1):
        start = cities[path[i]]
        end = cities[path[i + 1]]
        plt.plot([start[0], end[0]], [start[1], end[1]], 'r-', zorder=1)
    plt.title('TSP Visualization')
    plt.axis('equal')
    plt.show()