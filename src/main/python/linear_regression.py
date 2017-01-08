from __future__ import print_function
import sys
import time
import numpy as np
from pyspark import SparkContext

if __name__ == "__main__":

    sc = SparkContext(appName="PythonLR")
    D = 10
    p = 2
    iterations = 20
    N = 10
    if len(sys.argv)>1:
        N = int(sys.argv[1])
    if len(sys.argv)>2:
        iterations = int(sys.argv[2])
    if len(sys.argv)>3:
        D = int(sys.argv[3])
    print("N %d D %d iterations %d" %(N,D,iterations))

    points = sc.parallelize(range(1,N)).mapPartitions(lambda r: [np.random.ranf(size=(len(list(r)),D+p))])
    a = points.cache().first()
    start = time.time()
    alphaN = 0.01/N
    w = np.zeros(shape=(3,2))
    print("Initial w: " + str(w))

    def gradient(matrix, w):
        Y = matrix[:, 0:1]
        X = matrix[:, 2:]
        return alphaN * X.T.dot(X.dot(w)-Y)

    def add(x, y):
        x += y
        return x

    for i in range(iterations):
        print("On iteration %i" % (i + 1))
        w -= points.map(lambda m: gradient(m, w)).reduce(add)

    print("linear regression exec time %f" % (time.time()-start))
    print("Final w: " + str(w))

    sc.stop()
