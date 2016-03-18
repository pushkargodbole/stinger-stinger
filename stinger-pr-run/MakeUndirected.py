from collections import OrderedDict
from GraphIO import ReadGraph, WriteGraph
import sys

def MakeUndirected(DG):
    UG = OrderedDict()
    for u in DG:
        UG[u] = OrderedDict()
        for v in DG[u]:
            if v in UG and u in UG[v]:
                continue
                
            UG[u][v] = None

    return UG
    
if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise RuntimeError("Need filename of input graph")
    infile = sys.argv[1]
    outfile = infile.split('.')[0] + "_ug." + infile.split('.')[1]
    if len(sys.argv) > 2:
        graphtype = sys.argv[2]
    else:
        graphtype = "al"
        
    if len(sys.argv) > 3:
        header = int(sys.argv[3])
    else:
        header = 0
        
    if len(sys.argv) > 4:
        offset = int(sys.argv[4])
    else:
        offset = 0
        
    DG = ReadGraph(infile, graphtype, header, offset)
    print len(DG)
    UG = MakeUndirected(DG)
    WriteGraph(UG, outfile, graphtype)
