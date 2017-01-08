from __future__ import print_function

import sys
import time
import os

import numpy as np
from pyspark import SparkContext


def parseVector(line):
    return np.array([float(x) for x in line.split(',')])


def closestPoint(p, centers):
    bestIndex = 0
    closest = float("+inf")
    for i in range(len(centers)):
        tempDist = np.sum((p - centers[i]) ** 2)
        if tempDist < closest:
            closest = tempDist
            bestIndex = i
    return bestIndex


if __name__ == "__main__":

    sc = SparkContext(appName="PythonKMeans")
    file_name=os.environ['SCRATCH']+"/benchmark_data/kmeans_large.csv"
    lines = sc.textFile(file_name)
    data = lines.map(parseVector).cache()
    a = data.first()
    start = time.time()
    K = 5
    iterations = 10

    kPoints = data.takeSample(False, K, 1)

    for i in range(iterations):
        print("On iteration %i" % (i + 1))
        closest = data.map(
            lambda p: (closestPoint(p, kPoints), (p, 1)))
        pointStats = closest.reduceByKey(
            lambda p1_c1, p2_c2: (p1_c1[0] + p2_c2[0], p1_c1[1] + p2_c2[1]))
        newPoints = pointStats.map(
            lambda st: (st[0], st[1][0] / st[1][1])).collect()

        for (iK, p) in newPoints:
            kPoints[iK] = p

    print("kmeans exec time %f" % (time.time()-start))
    print("Final centers: " + str(kPoints))

    sc.stop()
