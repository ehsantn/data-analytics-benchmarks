from __future__ import print_function
import sys
import time
from random import random
from operator import add

from pyspark import SparkContext


if __name__ == "__main__":
    sc = SparkContext(appName="PythonPi")
    partitions = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    n = 100000 * partitions

    def f(_):
        x = random() * 2 - 1
        y = random() * 2 - 1
        return 1 if x ** 2 + y ** 2 < 1 else 0

    start = time.time()
    print("iterations: %d" % n)
    count = sc.parallelize(range(1, n + 1), partitions).map(f).reduce(add)
    print("exec time %f" % (time.time()-start))
    print("Pi is roughly %f" % (4.0 * count / n))

    sc.stop()
