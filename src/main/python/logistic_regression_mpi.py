import sys
import time
import numpy as np
import os
from mpi4py import MPI

if __name__ == "__main__":

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    pes = comm.Get_size()
    D = 10  # Number of dimensions
    iterations = 20
    N = 10
    if len(sys.argv)>1:
        N = int(sys.argv[1])
    if len(sys.argv)>2:
        iterations = int(sys.argv[2])
    if len(sys.argv)>3:
        D = int(sys.argv[3])
    chunk = int(N/pes)
    if rank==pes-1:
        chunk = N-(pes-1)*chunk
    if rank==0:
        print("N %d D %d iterations %d pes %d chunk %d" %(N,D,iterations,pes,chunk))

    # Initialize w to a random value
    X = np.random.ranf(size=(chunk,D))
    Y = np.random.ranf(size=(chunk))

    start = time.time()

    w = 2 * np.random.ranf(size=D) - 1
    if rank==0:
        print("Initial w: " + str(w))
    g_all = np.zeros(D)

    for i in range(iterations):
        g = ((1.0 / (1.0 + np.exp(-Y * X.dot(w))) - 1.0) * Y * X.T).sum(1)
        comm.Allreduce([g,MPI.DOUBLE],[g_all,MPI.DOUBLE])
        w -= g_all
    
    if rank==0:
        print("Final w: " + str(w))
        print("lr exec time %f" % (time.time()-start))
