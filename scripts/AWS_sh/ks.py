import sys
import time
import numpy as np
import os
from pyspark import SparkContext

n = 256000000
N = 3
b = 0.5
points = np.array([-1.0, 2.0, 5.0])

print("size n %d" %n)

X = sc.parallelize(range(1,n)).map(lambda r: np.random.random_sample())

X.cache().first()
start = time.time()

# Compute logistic regression gradient for a matrix of data points
def score(p):
	d = (-(p-points)**2)/(2*b**2)
	m = min(d)
	exps = m-np.log(b*N)+np.log(sum(np.exp(d-m)))
	return exps

def add(x, y):
	x += y
	return x

s = X.map(lambda m: score(m)).reduce(add)
print("Final s: " + str(s))

print("ks exec time %f" % (time.time()-start))


