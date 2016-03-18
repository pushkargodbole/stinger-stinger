# Randomize the ordering in an Edge List file
import sys
from GraphIO import ReadGraph, WriteGraph
import random

if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise RuntimeError("Need input adjacency-list file as argument")
    infile = sys.argv[1]
    outfile = infile.split('.')[0] + "_rand." + infile.split('.')[1]
    
    if len(sys.argv) > 2:
        header = int(sys.argv[2])
    else:
        header = 0
    
    if len(sys.argv) > 3:
        offset = int(sys.argv[3])
    else:
        offset = 0
    
    EL = ReadGraph(infile, "el", header, offset)
    el_list = []
    for u in EL:
        for v in EL[u]:
            el_list.append((u, v))
            
    random.shuffle(el_list)
    f = open(outfile, 'w')
    for elem in el_list:
        f.write(str(elem[0]) + ' ' + str(elem[1]) + '\n')
        
    f.close()
        

