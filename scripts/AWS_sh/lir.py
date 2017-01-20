from __future__ import print_function
import sys
import time
import numpy as np
from pyspark import SparkContext

D = 10
p = 4
iterations = 20
N = 256000000 
points = sc.parallelize(range(1,N)).mapPartitions(lambda r: [np.random.ranf(size=(len(list(r)),D+p))]).cache()
alphaN = 0.01/N
w = np.zeros(shape=(D,p))
print("Initial w: " + str(w))

def gradient(matrix, w):
	Y = matrix[:, 0:p]
	X = matrix[:, p:]
	return alphaN * X.T.dot(X.dot(w)-Y)

def add(x, y):
	 x += y
	 return x

w -= points.map(lambda m: gradient(m, w)).reduce(add)
start = time.time()

for i in range(iterations):
#print("On iteration %i" % (i + 1))
	w -= points.map(lambda m: gradient(m, w)).reduce(add)

print("linear regression exec time %f" % (time.time()-start))
print("Final w: " + str(w))

