from __future__ import print_function
import sys
import time
from pyspark.sql import SparkSession
from pyspark.sql import Row
from pyspark.sql.types import *



if __name__=="__main__":
    spark = SparkSession \
            .builder \
            .appName("UDF benchmark v2") \
            .config("spark.sql.autoBroadcastJoinThreshold", "-1") \
            .getOrCreate()
    sc = spark.sparkContext
    file_name = sys.argv[1]
    lines = sc.textFile(file_name)
    parts = lines.map(lambda l: l.split(","))
    data = parts.map(lambda p: (int(p[0]), float(p[1]), float(p[2])))
    schema = StructType([StructField("id",LongType(),True),\
              StructField("x",DoubleType(),True),\
              StructField("y",DoubleType(),True)])
    df = spark.createDataFrame(data, schema)
    df.cache().first()
    t1 = time.time()
    spark.udf.register("myudf", lambda x:2*x, DoubleType())
    df.registerTempTable("points")
    df1 = spark.sql("""SELECT
                    id AS id,
                    SUM(myudf(x))    AS sx,
                    SUM(myudf(y))    AS sy
                    FROM points 
                    GROUP BY id 
                    """)
    df1.cache().first()
    t2 = time.time()
    print("UDF python test v2 execution time %f" % (t2-t1))
    spark.stop()
