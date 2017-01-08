from __future__ import print_function
import sys
import time
import numpy as np
from pyspark import SparkContext

D = 10  # Number of dimensions

def readPointBatch(iterator):
    strs = list(iterator)
    matrix = np.zeros((len(strs), D + 1))
    for i, s in enumerate(strs):
        matrix[i] = np.fromstring(s.replace(',', ' '), dtype=np.float32, sep=' ')
    return [matrix]

if __name__ == "__main__":

    if len(sys.argv) != 3:
        print("Usage: logistic_regression <file> <iterations>", file=sys.stderr)
        exit(-1)

    sc = SparkContext(appName="PythonLR")
    points = sc.textFile(sys.argv[1]).mapPartitions(readPointBatch).cache()
    a = points.first()
    iterations = int(sys.argv[2])
    start = time.time()

    # Initialize w to a random value
    w = 2 * np.random.ranf(size=D).astype(np.float32) - 1
    print("Initial w: " + str(w))

    def gradient(matrix, w):
        Y = matrix[:, 0]    # point labels (first column of input file)
        X = matrix[:, 1:]   # point coordinates
        return ((1.0 / (1.0 + np.exp(-Y * X.dot(w))) - 1.0) * Y * X.T).sum(1)

    def add(x, y):
        x += y
        return x

    for i in range(iterations):
        print("On iteration %i" % (i + 1))
        w -= points.map(lambda m: gradient(m, w)).reduce(add)

    print("logistic regression exec time %f" % (time.time()-start))
    print("Final w: " + str(w))

    sc.stop()
