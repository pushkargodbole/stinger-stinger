import numpy as np
from block_matrix import gen_bdm, roc, gen_alist_file

def shift(i, step, n):
    if i-step >= 0:
        return i-step
    else:
        return n+i-step

def shiftup_bdm(matrix, step):
    n = matrix.shape[0]
    shifted_matrix = np.zeros(matrix.shape, dtype=int)
    for i in xrange(matrix.shape[0]):
        for j in xrange(matrix.shape[1]):
            shifted_matrix[shift(i, step, n)][shift(j, step, n)] = matrix[i][j]
    
    return shifted_matrix

def gen_add_del_matrices(old_matrix, new_matrix):
    add_matrix = np.where(new_matrix-old_matrix==1, 1, 0)
    del_matrix = np.where(new_matrix-old_matrix==-1, 1, 0)
    
    return add_matrix, del_matrix

if __name__ == "__main__":
    num_blocks = 10
    block_size = 20
    step = 10
    a = gen_bdm(num_blocks, block_size)
    roc(a, num_blocks, block_size)
    #a += np.identity(num_blocks*block_size, dtype=int)
    shifted_a = shiftup_bdm(a, step)
    add_a, del_a = gen_add_del_matrices(a, shifted_a)
    gen_alist_file(a, "roc_"+str(num_blocks)+"_"+str(block_size)+".csv")
    gen_alist_file(add_a, "roc_"+str(num_blocks)+"_"+str(block_size)+"_"+str(step)+"_shift_add.csv")
    gen_alist_file(del_a, "roc_"+str(num_blocks)+"_"+str(block_size)+"_"+str(step)+"_shift_del.csv")
    AL = ReadGraph("roc_"+str(num_blocks)+"_"+str(block_size)+".csv", "al", 1, 0)
    UDAL = MakeUndirected(AL)
    WriteGraph(UDAL, outfile, "el")
    """
    print a
    print
    print shifted_a
    print
    print add_a
    print
    print del_a
    """
