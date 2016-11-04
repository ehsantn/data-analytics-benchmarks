from __future__ import print_function
import sys 
import time

import numpy as np
import pandas as pd

if __name__=="__main__":

    file_name = sys.argv[1]
    df1 = pd.read_csv(file_name, header=None, names=['id','x','y'])
    t1 = time.time()
    df2 = df1[df1['id']<100]
    t2 = time.time()
    print("filter execution time %f" % (t2-t1))
