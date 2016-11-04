from __future__ import print_function
import sys
import time
from pyspark.sql import SparkSession
from pyspark.sql import Row
from pyspark.sql.types import *



if __name__=="__main__":
    spark = SparkSession \
            .builder \
            .appName("join benchmark") \
            .config("spark.sql.autoBroadcastJoinThreshold", "-1") \
            .getOrCreate()
    sc = spark.sparkContext
    file_name = sys.argv[1]
    lines = sc.textFile(file_name)
    parts = lines.map(lambda l: l.split(","))
    data1 = parts.map(lambda p: (int(p[0]), float(p[1]), float(p[2])))
    data2 = parts.map(lambda p: (int(p[0]), float(p[1]), float(p[2])))
    schema1 = StructType([StructField("id",LongType(),True),\
              StructField("x",DoubleType(),True),\
              StructField("y",DoubleType(),True)])
    schema2 = StructType([StructField("id",LongType(),True),\
              StructField("x1",DoubleType(),True),\
              StructField("y1",DoubleType(),True)])
    df1 = spark.createDataFrame(data1, schema1)
    df2 = spark.createDataFrame(data2, schema2)
    df1.cache().first()
    df2.cache().first()
    t1 = time.time()
    df3 = df1.join(df2, df1.id==df2.id, 'inner')
    df3.cache().first()
    t2 = time.time()
    print("join python test execution time %f" % (t2-t1))
    spark.stop()
