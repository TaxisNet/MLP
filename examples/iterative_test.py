import numpy as np
import networkx as nx
# add to PYTHONPATH
import sys
sys.path.append('/home/taxis/Documents/MLP/build')
import mlp_py


def iterative_solutions(soloutions):
    for i in range(len(soloutions)-1):
        if solutions[i][1:] != solutions[i+1]:
            print(f"Solution {i} and Solution {i+1} are different.")
            return False
    print("All solutions are the same.")
    return True


def get_distance_matrix(graph, start_node=None):
    n = graph.number_of_nodes()
    
    if start_node is not None:
        # Reorder the graph so that the start_node is at index 0
        nodes = list(graph.nodes())
        if start_node in nodes:
            nodes.remove(start_node)
            nodes.insert(0, start_node)
        else:
            raise ValueError(f"Start node {start_node} not found in the graph.")
    
    mapping = {idx: node for idx, node in enumerate(nodes)}
    distance_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            pos_i = graph.nodes[mapping[i]]['position']
            pos_j = graph.nodes[mapping[j]]['position']
            dist = np.linalg.norm(pos_i - pos_j)  # Euclidean distance between nodes
            distance_matrix[i][j] = dist
            distance_matrix[j][i] = dist  # Symmetric matrix
    return distance_matrix, mapping

def get_weights(graph, mapping):
    n = graph.number_of_nodes()
    weights = np.zeros(n)
    for  idx, node in mapping.items():
        weights[idx] = graph.nodes[node]['weight']
    weights /= np.sum(weights)  # Normalize weights to sum to 1
    return weights


np.random.seed(100)  # For reproducibility
n = 50
# Make a graph with n nodes and random weights
weights = np.random.rand(n)  # Random weights for each node
weights /= np.sum(weights)  # Normalize weights to sum to 1
positions = np.random.rand(n, 2) * 100  # Random positions for each node in a 100x100 space
G = nx.Graph()
for i in range(n):
    G.add_node(i, weight=weights[i], position = positions[i])  # Random weight for each node


distance_matrix, mapping = get_distance_matrix(G, start_node=0)
weights = get_weights(G, mapping)


solutions = []
costs = []

for i in range(n):
    print(f"Iteration {i+1}/{n}")
    cost, solution = mlp_py.run(distance_matrix, weights, verbose=False)
    solutions.append([mapping[node] for node in solution])  # Store the solution with original node labels
    costs.append(cost)
    mapped_solution = [mapping[node] for node in solution]
    G.remove_node(mapped_solution[0])  # Remove the node from the graph
    if i == n - 1:
        mapped_solution = [mapping[node] for node in solution[1:]]  # Adjust for the last iteration
        break  # No need to update for the last iteration
    start = mapped_solution[1]
    distance_matrix, mapping = get_distance_matrix(G, start_node=start)  # Update the distance matrix
    weights = get_weights(G, mapping)  # Update the weights


print(iterative_solutions(solutions))