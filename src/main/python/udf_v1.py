from __future__ import print_function
import sys
import time
from pyspark.sql import SparkSession
from pyspark.sql import Row
from pyspark.sql.types import *



if __name__=="__main__":
    spark = SparkSession \
            .builder \
            .appName("UDF benchmark v1") \
            .config("spark.sql.autoBroadcastJoinThreshold", "-1") \
            .getOrCreate()
    sc = spark.sparkContext
    file_name = sys.argv[1]
    lines = sc.textFile(file_name)
    parts = lines.map(lambda l: l.split(","))
    data = parts.map(lambda p: (p[0], p[1], p[2]))
    schema = StructType(StructField("id",LongType(),True),\
              StructField("x",DoubleType(),True),\
              StructField("y",DoubleType(),True))
    df = spark.createDataFrame(data, schema)
    df.cache().first()
    t1 = time.time()
    df.registerTempTable("points")
    df1 = spark.sql("""SELECT
                    id AS id,
                    SUM(2*x)    AS sx,
                    SUM(2*y)    AS sy
                    FROM points 
                    GROUP BY id 
                    """)
    df1.cache().first()
    t2 = time.time()
    print("UDF test v1 execution time %f" % (t2-t1))
    spark.stop()
