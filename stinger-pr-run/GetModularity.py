# Given edge list and community labeling for each node, compute modularity

import sys

def GetModularity(edge_list, labels, offset=0):
    graph = {}
    for node in xrange(len(labels)):
        graph[node] = set()
        
    for edge in edge_list:
        graph[(edge[0]-offset)].add((edge[1]-offset))
        graph[(edge[1]-offset)].add((edge[0]-offset))
    
    mod = 0.
    for node in xrange(len(labels)):
        for nbr in xrange(len(labels)):
            if labels[(node+offset)] == labels[(nbr+offset)]:
                mod += (1 if nbr in graph[node] else 0) - float(len(graph[node])*len(graph[nbr]))/(2*len(edge_list))
                #print node, nbr, (1 if nbr in graph[node] else 0), len(graph[node]), len(graph[nbr]), (1 if nbr in graph[node] else 0) - float(len(graph[node])*len(graph[nbr]))/(2*len(edge_list))
                
    mod /= 2*len(edge_list)
    
    return mod
    
if __name__ == "__main__":
    if len(sys.argv) < 3:
        raise RuntimeError("Requires at least 2 arguments <edge_list_filename> <community_labels_filename> [<node_offset>]")
    
    if len(sys.argv) == 4:
        offset = int(sys.argv[3])
    else:
        offset = 0
        
    edge_list = []
    with open(sys.argv[1], 'r') as el:
        for line in el:
            edge = (int(line.strip().split()[0]), int(line.strip().split()[1]))
            edge_list.append(edge)
    
    labels = {}        
    with open(sys.argv[2], 'r') as l:
        for line in l:
            node = int(line.strip().split()[0])
            label = int(line.strip().split()[1])
            labels[node] = label
            
    print GetModularity(edge_list, labels, offset)
                    
                    
