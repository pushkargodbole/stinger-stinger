# Convert graph file from Edge List to Adjacency List
import sys
from GraphIO import ReadGraph, WriteGraph

if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise RuntimeError("Need input adjacency-list file as argument")
    infile = sys.argv[1]
    outfile = infile.split('.')[0].replace("_el", "") + "_al." + infile.split('.')[1]
    
    if len(sys.argv) > 2:
        header = int(sys.argv[2])
    else:
        header = 0
    
    EL = ReadGraph(infile, "el", header)
    WriteGraph(EL, outfile, "al")
    
