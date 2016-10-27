
import org.apache.spark.{SparkConf, SparkContext}
import org.apache.spark.sql.DataFrame
import org.apache.spark.sql.Row
import org.apache.spark.sql.types._
import org.apache.spark.sql.SQLContext
import org.apache.spark.sql.functions._
import scala.language.existentials

import org.apache.spark.sql.catalyst.analysis.UnresolvedRelation
import org.apache.spark.sql.catalyst.TableIdentifier
import org.apache.spark.sql.execution.joins._
import org.apache.spark.ml.linalg.{Vector, Vectors}
import org.apache.spark.ml.feature.VectorAssembler

import java.lang.management.ManagementFactory
import scala.collection.JavaConversions._

/**
  *  q05
command to run
./bin/spark-submit --num-executors 4  --jars /home/whassan/Downloads/commons-csv-1.2.jar,/home/whassan/Downloads/spark-csv_2.10-1.4.0.jar  --class Query05  ~/hps/query-examples/target/scala-2.10/query05_2.10-0.1.jar "/home/whassan/tmp/csv/store_sales_sanitized_200_lines.csv"  "/home/whassan/tmp/csv/item_sanitized_200_lines.csv"

  store_sales table absolute path = args[0]
  item table absolute path = args[1]

  */
object Query05 {
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
    val query_05_sql_path = System.getProperty("user.home") + "/pse-hpc/spark-sql-query-tests/src/main/scala/q05.sql"

    val table_web_clickstreams_path = args(0)
    val table_item_path = args(1)
    val table_customer_path = args(2)
    val table_customer_demographics_path = args(3)

    val schema_web_clickstreams = StructType(Array(
      StructField("wcs_item_sk",IntegerType,true),
      StructField("wcs_user_sk",IntegerType,true)))

    val df_web_clickstreams = sqlContext.read.schema(schema_web_clickstreams).csv(table_web_clickstreams_path)
    df_web_clickstreams.cache().count()

    val schema_item = StructType(Array(
      StructField("i_item_sk",IntegerType,true),
      StructField("i_category_id",IntegerType,true),
      StructField("i_category",LongType,true)))
    val df_item = sqlContext.read.schema(schema_item).csv(table_item_path)
    df_item.cache().first()

    val schema_customer = StructType(Array(
      StructField("c_customer_sk",IntegerType,true),
      StructField("c_current_cdemo_sk",IntegerType,true)))
    val df_customer = sqlContext.read.schema(schema_customer).csv(table_customer_path)
    df_customer.cache().first()

    val schema_customer_demographics = StructType(Array(
      StructField("cd_demo_sk",IntegerType,true),
      StructField("cd_gender",IntegerType,true),
      StructField("cd_education_status",IntegerType,true)))
    val df_customer_demographics = sqlContext.read.schema(schema_customer_demographics).csv(table_customer_demographics_path)
    df_customer_demographics.cache().first()

    df_web_clickstreams.registerTempTable("web_clickstreams")
    df_item.registerTempTable("item")
    df_customer.registerTempTable("customer")
    df_customer_demographics.registerTempTable("customer_demographics")

    val t0 = System.currentTimeMillis
    val lines = scala.io.Source.fromFile(query_05_sql_path).mkString
    val fin  = sqlContext.sql(lines)
    // fin.cache().first()
    val assembler = new VectorAssembler()
      .setInputCols(Array("clicks_in_category","college_education","male","clicks_in_1","clicks_in_2","clicks_in_3","clicks_in_4","clicks_in_5","clicks_in_6","clicks_in_7"))
      .setOutputCol("features")
    val ds = assembler.transform(fin)
    ds.cache.first
    val t1 = System.currentTimeMillis
    println("****** Query 05 time(ms) took: " + (t1 - t0).toFloat/1000.0)
    fin.show()
    printExecutionPlan(fin)
    println(":Done with Query 05")
  }
}
