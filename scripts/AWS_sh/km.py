from __future__ import print_function
import sys
import time
import os
import numpy as np
from pyspark import SparkContext

def closestPoint(p, centers):
    bestIndex = 0
    closest = float("+inf")
    for i in range(len(centers)):
        tempDist = np.sum((p - centers[i]) ** 2)
        if tempDist < closest:
            closest = tempDist
            bestIndex = i
    return bestIndex



D = 10
N = 64000000
K = 5
iterations = 20
print("N %d D %d iterations %d K %d" %(N,D,iterations,K))

data = sc.parallelize(range(1,N)).map(lambda r: np.random.ranf(size=D))
a = data.cache().first()
start = time.time()

kPoints = data.takeSample(False, K, 1)

for i in range(iterations):
	#print("On iteration %i" % (i + 1))
	closest = data.map(lambda p: (closestPoint(p, kPoints), (p, 1)))
	pointStats = closest.reduceByKey(lambda p1_c1, p2_c2: (p1_c1[0] + p2_c2[0], p1_c1[1] + p2_c2[1]))
	newPoints = pointStats.map(lambda st: (st[0], st[1][0] / st[1][1])).collect()
	for (iK, p) in newPoints:
		kPoints[iK] = p

print("kmeans exec time %f" % (time.time()-start))
print("Final centers: " + str(kPoints))

