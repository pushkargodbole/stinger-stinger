import matplotlib.pyplot as plt
import matplotlib
import sys

def Plot(filename, m1, m2=None, m3=None):
    f = open(filename, 'r')
    i = 0
    header = None
    c1 = None
    c2 = None
    c3 = None
    C1 = []
    C2 = []
    C3 = []
    for line in f:
        if i == 0:
            header = line.strip().split()
            for idx, field in enumerate(header):
                if field == m1:
                    c1 = idx
                if m2!=None and field == m2:
                    c2 = idx
                if m3!=None and field == m3:
                    c3 = idx
        else:
            elems = line.strip().split()
            C1.append(float(elems[c1]))
            if c2 != None:
                C2.append(float(elems[c2]))
            if c3 != None:
                C3.append(float(elems[c3]))
        i+=1
    
    matplotlib.rcParams.update({'font.size': 16.5})    
    plt.plot(C1, 'ro')
    p1, = plt.plot(C1, 'r-')
    if c2 != None:
        plt.plot(C2, 'bo')
        p2, = plt.plot(C2, 'b-')
    if c3 != None:
        plt.plot(C3, 'go')
        p3, = plt.plot(C3, 'g-')
    plt.xlabel('Batch no.')
    plt.ylabel(header[c1].split('(')[0])
    if c2 == None:
        plt.title(filename + '\n\n' + header[c1])
        plt.legend([p1], [header[c1]])
    elif c3 == None:
        plt.title(filename + '\n\n' + header[c1] + ' vs ' + header[c2])
        plt.legend([p1, p2], [header[c1], header[c2]], loc=1)
    else:
        plt.title(filename + '\n\n' + header[c1] + ' vs ' + header[c2] + ' vs ' + header[c3])
        plt.legend([p1, p2, p3], [header[c1], header[c2], header[c3]])
    plt.show()
    
if __name__ == "__main__":
    filename = sys.argv[1]
    m1 = sys.argv[2]
    m2 = None
    m3 = None
    if len(sys.argv) > 3:
        m2 = sys.argv[3]
    if len(sys.argv) > 4:
        m3 = sys.argv[4]
    
    Plot(filename, m1, m2, m3)
