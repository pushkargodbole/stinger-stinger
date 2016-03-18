# Convert graph file from Adjacency List to Edge List
import sys
from GraphIO import ReadGraph, WriteGraph

if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise RuntimeError("Need input adjacency-list file as argument")
    infile = sys.argv[1]
    outfile = infile.split('.')[0].replace("_al", "") + "_el." + infile.split('.')[1]
    
    if len(sys.argv) > 2:
        header = int(sys.argv[2])
    else:
        header = 0
    
    if len(sys.argv) > 3:
        offset = int(sys.argv[3])
    else:
        offset = 0
    
    AL = ReadGraph(infile, "al", header, offset)
    WriteGraph(AL, outfile, "el")
        
