import org.apache.spark.{SparkConf, SparkContext}
import org.apache.spark.SparkContext._
import org.apache.spark.sql.SparkSession
import org.apache.spark.sql.DataFrame
import org.apache.spark.sql._
import org.apache.spark.sql.types._
import org.apache.spark.sql.SQLContext
import org.apache.spark.sql.functions._
import scala.language.existentials

import org.apache.spark.sql.catalyst.analysis.UnresolvedRelation
import org.apache.spark.sql.catalyst.TableIdentifier
import org.apache.spark.sql.execution.joins._

import org.apache.spark.ml.clustering.{KMeansModel, KMeans}
import org.apache.spark.ml.linalg.{Vector, Vectors}
import org.apache.spark.ml.feature.VectorAssembler
import org.apache.spark.rdd.RDD
import org.apache.spark.sql._
import org.apache.spark.sql.types._
import org.apache.spark.sql.catalyst.encoders.ExpressionEncoder
// import spark.implicits._
import scala.language.reflectiveCalls
import java.lang.management.ManagementFactory
import scala.collection.JavaConversions._
/**
 *  q26
 command to run
 ./bin/spark-submit --num-executors 4  --jars /home/whassan/Downloads/commons-csv-1.2.jar,/home/whassan/Downloads/spark-csv_2.10-1.4.0.jar  --class Query26  ~/hps/query-examples/target/scala-2.10/query26_2.10-0.1.jar "/home/whassan/tmp/csv/store_sales_sanitized_200_lines.csv"  "/home/whassan/tmp/csv/item_sanitized_200_lines.csv"

store_sales table absolute path = args[0]
item table absolute path = args[1]

*/

object Query25 {
  // print the execution plan
  def printExecutionPlan(fin: DataFrame){
    println(fin.queryExecution.logical.numberedTreeString)
    println("\n===================================\n")
    println(fin.queryExecution.optimizedPlan.numberedTreeString)
    println("\nExecuted Plan=====================\n")
    println(fin.queryExecution.executedPlan.numberedTreeString)
    println("\nSpark Plan=====================\n")
    println(fin.queryExecution.sparkPlan.numberedTreeString)
    println("\nStatistics=====================\n")
    println(fin.queryExecution.analyzed.statistics.sizeInBytes)
    println(fin.queryExecution.toString)
  }
  def main(args: Array[String]) {
    val spark = SparkSession
      .builder()
      .appName("Q25")
      .config("spark.sql.autoBroadcastJoinThreshold", "-1")
      .getOrCreate()

    import spark.implicits._
    val table_store_sales_path = args(0)
    val table_web_sales_path = args(1)
    // Change this variable where sql is

    val schema_store_sales = StructType(Array(
      StructField("ss_sold_date_sk", IntegerType,true),
      StructField("ss_customer_sk", IntegerType,true),
      StructField("ss_ticket_number", IntegerType,true),
      StructField("ss_net_paid", DoubleType,true)))
    val df_store_sales = spark.read.schema(schema_store_sales).csv(table_store_sales_path)

    val schema_web_sales = StructType(Array(
      StructField("ws_sold_date_sk", IntegerType,true),
      StructField("ws_bill_customer_sk", IntegerType,true),
      StructField("ws_order_number", IntegerType,true),
      StructField("ws_net_paid", DoubleType,true)))
    val df_web_sales = spark.read.schema(schema_web_sales).csv(table_web_sales_path)

    df_store_sales.registerTempTable("store_sales_table")
    df_web_sales.registerTempTable("web_sales_table")
    // collect() fails so using first()
    df_store_sales.cache().first()
    df_web_sales.cache().first()
    // Starting time
    val t0 = System.currentTimeMillis
    val df1 = spark.sql("""SELECT
        ss_customer_sk          AS cid,
        count(distinct ss_ticket_number)   AS frequency,
        max(ss_sold_date_sk)               AS most_recent_date,
        SUM(ss_net_paid)                   AS amount
        FROM store_sales_table
        WHERE ss_sold_date_sk > 33000
        GROUP BY ss_customer_sk
        """ )
    val df2 = spark.sql("""SELECT
        ws_bill_customer_sk                AS cid,
        count(distinct ws_order_number)    AS frequency,
        max(ws_sold_date_sk)               AS most_recent_date,
        SUM(ws_net_paid)                   AS amount
        FROM web_sales_table
        WHERE ws_sold_date_sk > 33000
        GROUP BY ws_bill_customer_sk
        """ )
    val df3 = df1.union(df2)
    df3.registerTempTable("agg_sales")
    val df4 = spark.sql("""SELECT
        cid as cid,
        CASE WHEN 37621 - max(most_recent_date) < 60 THEN 1.0 ELSE 0.0 END 
        AS recency, -- 37621 == 2003-01-02
        SUM(frequency) AS frequency, --total frequency
        SUM(amount)    AS totalspend --total amount
        FROM agg_sales
        GROUP BY cid 
        """ )


    // head to intiate lazy evaluation
    //df4.cache.head
    // From spark website, there should be a good way
    val assembler = new VectorAssembler()
      .setInputCols(Array("cid","recency","frequency","totalspend"))
      .setOutputCol("features")
    val ds = assembler.transform(df4)
    ds.cache.first
    val t1 = System.currentTimeMillis
    // Clusters = 8  and Iterations 20
    val means = new KMeans().setK(8).setMaxIter(20)
    means.setInitMode("random").setSeed(675234312453645L)
    val clusterModel= means.fit(ds)
    // head is called to initiate the lazy evaluation
    clusterModel.clusterCenters.head
    // Measure time
    println("****** Query 25 time(s) took: " + (t1 - t0).toFloat / 1000)
    println(":Done with Query 25")
  }
}
