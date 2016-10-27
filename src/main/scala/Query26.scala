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

object Query26 {
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
    val spark = SparkSession
      .builder()
      .appName("Q26")
      .config("spark.sql.autoBroadcastJoinThreshold", "-1")
      .getOrCreate()

    import spark.implicits._
    val table_store_slaes_path = args(0)
    val table_item_path = args(1)
    // Change this variable where sql is
    val query_26_sql_path = System.getProperty("user.home") + "/pse-hpc/spark-sql-query-tests/src/main/scala/q26.sql"

    val schema_store_sales = StructType(Array(
      StructField("ss_item_sk",IntegerType,true),
      StructField("ss_customer_sk",IntegerType,true)))
    val df_store_sales = spark.read.schema(schema_store_sales).csv(table_store_slaes_path)

    val schema_item = StructType(Array(
      StructField("i_item_sk",IntegerType,true),
      StructField("i_class_id",IntegerType,true),
      StructField("i_category",LongType,true)))

    val df_item = spark.read.schema(schema_item).csv(table_item_path)
    df_store_sales.registerTempTable("store_sales_table")
    df_item.registerTempTable("item_table")
    // collect() fails so using first()
    df_store_sales.cache().first()
    df_item.cache().first()
    // Starting time
    val t0 = System.currentTimeMillis
    val lines = scala.io.Source.fromFile(query_26_sql_path).mkString
    val fin  = spark.sql(lines)
    // head to intiate lazy evaluation
    // fin.cache.head
    // From spark website, there should be a good way
    val assembler = new VectorAssembler()
      .setInputCols(Array("id1", "id2", "id3","id4","id5","id6","id7","id8","id9","id10","id11","id12","id13","id14","id15"))
      .setOutputCol("features")
    val ds = assembler.transform(fin)
    ds.cache.first
    val t1 = System.currentTimeMillis
    // Clusters = 8  and Iterations 20
    val means = new KMeans().setK(8).setMaxIter(20)
    means.setInitMode("random").setSeed(675234312453645L)
    val clusterModel= means.fit(ds)
    // head is called to initiate the lazy evaluation
    clusterModel.clusterCenters.head
    // Measure time
    println("****** Query 26 time(s) took: " + (t1 - t0).toFloat / 1000)
    println(":Done with Query 26")
  }
}
