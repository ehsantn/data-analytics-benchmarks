from __future__ import print_function
import sys 
import time

import numpy as np
import pandas as pd

if __name__=="__main__":

    file_name = sys.argv[1]
    df1 = pd.read_csv(file_name, header=None, names=['id','x','y'])
    df2 = pd.read_csv(file_name, header=None, names=['id','x1','y1'])
    t1 = time.time()
    df3 = pd.merge(df1, df2) 
    t2 = time.time()
    print("join execution time %f" % (t2-t1))
