import sys
import time
import numpy as np
import os
from pyspark import SparkContext


if __name__ == "__main__":

    sc = SparkContext(appName="LR")
    D = 10  # Number of dimensions
    iterations = 20
    N = 10
    if len(sys.argv)>1:
        N = int(sys.argv[1])
    if len(sys.argv)>2:
        iterations = int(sys.argv[2])
    if len(sys.argv)>3:
        D = int(sys.argv[3])
    print("N %d D %d iterations %d" %(N,D,iterations))

    points = sc.parallelize(range(1,N)).mapPartitions(lambda r: [np.random.ranf(size=(len(list(r)),D+1))])

    points.cache().first()
    start = time.time()

    w = 2 * np.random.ranf(size=D) - 1
    print("Initial w: " + str(w))

    # Compute logistic regression gradient for a matrix of data points
    def gradient(matrix, w):
        Y = matrix[:, 0]    # point labels (first column of input file)
        X = matrix[:, 1:]   # point coordinates
        # For each point (x, y), compute gradient function, then sum these up
        return ((1.0 / (1.0 + np.exp(-Y * X.dot(w))) - 1.0) * Y * X.T).sum(1)

    def add(x, y):
        x += y
        return x

    for i in range(iterations):
        #print("On iteration %i" % (i + 1))
        w -= points.map(lambda m: gradient(m, w)).reduce(add)

    print("Final w: " + str(w))

    print("lr exec time %f" % (time.time()-start))
    sc.stop()
