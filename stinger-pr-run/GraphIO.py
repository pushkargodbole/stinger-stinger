from collections import OrderedDict

def ReadGraph(filename, graphtype="al", header=0, offset=0):
    if graphtype == "al":
        gt = "adjacency-list"
    elif graphtype == "el":
        gt = "edge-list"
    else:
        raise ValueError("graphtype has to be 'al'/'el'")
        
    print "Reading the {0} graph file '{1}', with {2} header {4} and node offset of {3}".format(gt, filename, header, offset, "line" if (header==1) else "lines")
    f = open(filename, 'r')
    u = 0
    Graph = OrderedDict()
    if graphtype == "al":
        for line in f:
            if u == header:
                u = 0
                header = -1
            if header == -1:    
                elems = line.strip().split()
                Graph[u+offset] = OrderedDict()
                for elem in elems:
                    v = int(elem.strip())
                    Graph[u+offset][v] = None
            u+=1
    else:
        for line in f:
            if u == header:
                u = 0
                header = -1
            if header == -1:    
                elems = line.strip().split()
                if int(elems[0]) not in Graph:
                    Graph[int(elems[0])] = OrderedDict()
                Graph[int(elems[0])][int(elems[1])] = None
            u+=1
            
    f.close()
    
    return Graph
    
def WriteGraph(Graph, filename, graphtype="al"):
    if graphtype == "al":
        gt = "adjacency-list"
    elif graphtype == "el":
        gt = "edge-list"
    else:
        raise ValueError("graphtype has to be 'al'/'el'")
    
    offset = min(Graph.keys())
    print "Writing the {0} graph to file '{1}', with node offset of {2}".format(gt, filename, offset)
    f = open(filename, 'w')
    if graphtype == "al":
        for u in Graph:
            line = ''
            for v in Graph[u]:
                line += str(v) + ' '
            
            line += '\n';
            f.write(line)
    else:
        for u in Graph:
            for v in Graph[u]:
                line = str(u) + ' ' + str(v) + '\n'
                f.write(line)

    f.close()
