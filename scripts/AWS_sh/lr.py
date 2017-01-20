import sys
import time
import numpy as np
from pyspark import SparkContext

D = 10  # Number of dimensions
iterations = 20
N = 256000000
print("N %d D %d iterations %d" %(N,D,iterations))

# Initialize w to a random value
points = sc.parallelize(range(1,N)).mapPartitions(lambda r: [np.random.ranf(size=(len(list(r)),D+1))]).cache()
def gradient(matrix, w): 
	Y = matrix[:, 0]    # point labels (first column of input file)
	X = matrix[:, 1:]   # point coordinates
	return ((1.0 / (1.0 + np.exp(-Y * X.dot(w))) - 1.0) * Y * X.T).sum(1)

def add(x, y): 
	x += y
	return x

w = 2 * np.random.ranf(size=D) - 1

print("Initial w: " + str(w))

start = time.time()
w -= points.map(lambda m: gradient(m, w)).reduce(add) 

for i in range(iterations):
	print("On iteration %d" % (i + 1))
	w -= points.map(lambda m: gradient(m, w)).reduce(add) 

print("Final w: " + str(w))

print("lr exec time %f" % (time.time()-start))

