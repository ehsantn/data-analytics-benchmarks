import org.apache.spark.{SparkConf, SparkContext}
import org.apache.spark.SparkContext._
import org.apache.spark.sql.DataFrame
import org.apache.spark.sql.Row
import org.apache.spark.sql.types._
import org.apache.spark.sql.SQLContext
import org.apache.spark.sql.functions._
import scala.language.existentials

import org.apache.spark.sql.catalyst.analysis.UnresolvedRelation
import org.apache.spark.sql.catalyst.TableIdentifier
import org.apache.spark.sql.execution.joins._

import org.apache.spark.mllib.clustering.{KMeansModel, KMeans}
import org.apache.spark.mllib.linalg.{Vector, Vectors}
import org.apache.spark.mllib.regression.LabeledPoint
import org.apache.spark.rdd.RDD

import java.lang.management.ManagementFactory
import scala.collection.JavaConversions._
/**
 *  q26
 command to run
 ./bin/spark-submit --num-executors 4  --jars /home/whassan/Downloads/commons-csv-1.2.jar,/home/whassan/Downloads/spark-csv_2.10-1.4.0.jar  --class Query26  ~/hps/query-examples/target/scala-2.10/query26_2.10-0.1.jar "/home/whassan/tmp/csv/store_sales_sanitized_200_lines.csv"  "/home/whassan/tmp/csv/item_sanitized_200_lines.csv"

store_sales table absolute path = args[0]
item table absolute path = args[1]

*/

object Query26_mllib {
  // print the execution plance
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
    val sparkConf = new SparkConf().setAppName("Query")
    val sc = new SparkContext(sparkConf)
    val sqlContext = new org.apache.spark.sql.SQLContext(sc)
    sqlContext.setConf("spark.sql.autoBroadcastJoinThreshold","-1")

    val table_store_slaes_path = args(0)
    val table_item_path = args(1)
    // Change this variable where sql is
    val query_26_sql_path = "/home/whassan/spark-sql-query-tests/src/main/scala/q26.sql"

    val schema_store_sales = StructType(Array(
      StructField("ss_item_sk",IntegerType,true),
      StructField("ss_customer_sk",IntegerType,true)))
    val df_store_sales = sqlContext.read.schema(schema_store_sales).csv(table_store_slaes_path)

    val schema_item = StructType(Array(
      StructField("i_item_sk",IntegerType,true),
      StructField("i_class_id",IntegerType,true),
      StructField("i_category",StringType,true)))

    val df_item = sqlContext.read.schema(schema_item).csv(table_item_path)
    df_store_sales.registerTempTable("store_sales_table")
    df_item.registerTempTable("item_table")
    // collect() fails so using first()
    df_store_sales.cache().first()
    df_item.cache().first()
    // Starting time
    val t0 = System.currentTimeMillis
    val lines = scala.io.Source.fromFile(query_26_sql_path).mkString
    val fin  = sqlContext.sql(lines)
    // head to intiate lazy evaluation
    fin.cache.head
    val t1 = System.currentTimeMillis
    // From spark website, there should be a good way
    val vectors = fin.rdd.map(s => Vectors.dense((s.get(0).toString)toDouble,(s.get(1).toString)toDouble,  (s.get(2).toString)toDouble, (s.get(3).toString)toDouble,(s.get(4).toString)toDouble,(s.get(5).toString)toDouble,(s.get(6).toString)toDouble,(s.get(7).toString)toDouble,(s.get(8).toString)toDouble,(s.get(9).toString)toDouble,(s.get(10).toString)toDouble,(s.get(11).toString)toDouble,(s.get(12).toString)toDouble,(s.get(13).toString)toDouble,(s.get(14).toString)toDouble, (s.get(15).toString)toDouble)).cache
    // Clusters = 8  and Iterations 20
    val means = new KMeans().setK(8).setMaxIterations(20)
    means.setInitializationMode(KMeans.RANDOM).setSeed(675234312453645L)
    val clusterModel= means.run(vectors)
    // head is called to initiate the lazy evaluation
    clusterModel.clusterCenters.head
    // Measure time
    println("****** Query 26 time(s) took: " + (t1 - t0).toFloat / 1000)
    println(":Done with Query 26")
  }
}
