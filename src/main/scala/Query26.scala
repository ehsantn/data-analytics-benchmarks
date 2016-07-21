
import org.apache.spark.{SparkConf, SparkContext}
import org.apache.spark.sql.Row
import org.apache.spark.sql.types._
import org.apache.spark.sql.SQLContext
import org.apache.spark.sql.functions._
import com.databricks.spark.csv._
import scala.language.existentials

import org.apache.spark.sql.catalyst.analysis.UnresolvedRelation
import org.apache.spark.sql.catalyst.TableIdentifier
import org.apache.spark.sql.execution.joins._
/**
  *  q26
command to run
./bin/spark-submit --num-executors 4  --jars /home/whassan/Downloads/commons-csv-1.2.jar,/home/whassan/Downloads/spark-csv_2.10-1.4.0.jar  --class Query26  ~/hps/query-examples/target/scala-2.10/query26_2.10-0.1.jar "/home/whassan/tmp/csv/store_sales_sanitized_200_lines.csv"  "/home/whassan/tmp/csv/item_sanitized_200_lines.csv"

  store_sales table absolute path = args[0]
  item table absolute path = args[1]

  */
object Query26 {
  def main(args: Array[String]) {
    val sparkConf = new SparkConf().setAppName("Query")
    val sc = new SparkContext(sparkConf)
    val sqlContext = new org.apache.spark.sql.SQLContext(sc)
    //sqlContext.setConf("spark.sql.autoBroadcastJoinThreshold","-1")
    //sqlContext.setConf("spark.broadcast.factory","None")
    //sqlContext.setConf("spark.broadcast.blockSize","0")
    //sqlContext.setConf("spark.driver.maxResultSize","16g")
    val table_store_slaes_path = args(0)
    val table_item_path = args(1)
    // uncomment following for complete dataset
    // val schema_store_sales = StructType(Array(
    //   StructField("ss_sold_date_sk", IntegerType,true),
    //   StructField("ss_sold_time_sk",IntegerType,true),
    //   StructField("ss_item_sk",IntegerType,true),
    //   StructField("ss_customer_sk",IntegerType,true),
    //   StructField("ss_cdemo_sk",IntegerType,true),
    //   StructField("ss_hdemo_sk",IntegerType,true),
    //   StructField("ss_addr_sk",IntegerType,true),
    //   StructField("ss_store_sk",IntegerType,true),
    //   StructField("ss_promo_sk",IntegerType,true),
    //   StructField("ss_ticket_number",IntegerType,true),
    //   StructField("ss_quantity",IntegerType,true),
    //   StructField("ss_wholesale_cost",DoubleType,true),
    //   StructField("ss_list_price",DoubleType,true),
    //   StructField("ss_sales_price",DoubleType,true),
    //   StructField("ss_ext_discount_amt",DoubleType,true),
    //   StructField("ss_ext_sales_price",DoubleType,true),
    //   StructField("ss_ext_wholesale_cost",DoubleType,true),
    //   StructField("ss_ext_list_price",DoubleType,true),
    //   StructField("ss_ext_tax",DoubleType,true),
    //   StructField("ss_coupon_amt",DoubleType,true),
    //   StructField("ss_net_paid",DoubleType,true),
    //   StructField("ss_net_paid_inc_tax",DoubleType,true),
    //   StructField("ss_net_profit",DoubleType,true)))


    val schema_store_sales = StructType(Array(
      StructField("ss_item_sk",IntegerType,true),
      StructField("ss_customer_sk",IntegerType,true)))
    val df_store_sales = sqlContext.read.format("com.databricks.spark.csv").option("header", "true").schema(schema_store_sales).load(table_store_slaes_path)

    // uncomment following for complete dataset
    // val schema_item = StructType(Array(
    //   StructField("i_item_sk",IntegerType,true),
    //   StructField("i_item_id",StringType,true),
    //   StructField("i_rec_start_date",StringType,true),
    //   StructField("i_rec_end_date",StringType,true),
    //   StructField("i_item_desc",StringType,true),
    //   StructField("i_current_price",DoubleType,true),
    //   StructField("i_wholesale_cost",DoubleType,true),
    //   StructField("i_brand_id",IntegerType,true),
    //   StructField("i_brand",StringType,true),
    //   StructField("i_class_id",IntegerType,true),
    //   StructField("i_class",StringType,true),
    //   StructField("i_category_id",IntegerType,true),
    //   StructField("i_category",StringType,true),
    //   StructField("i_manufact_id",IntegerType,true),
    //   StructField("i_manufact",StringType,true),
    //   StructField("i_size",StringType,true),
    //   StructField("i_formulation",StringType,true),
    //   StructField("i_color",StringType,true),
    //   StructField("i_units",StringType,true),
    //   StructField("i_container",StringType,true),
    //   StructField("i_manager_id",IntegerType,true),
    //   StructField("i_product_name",StringType,true)))

    val schema_item = StructType(Array(
      StructField("i_item_sk",IntegerType,true),
      StructField("i_class_id",IntegerType,true),
      StructField("i_category",StringType,true)))

    val df_item = sqlContext.read.format("com.databricks.spark.csv").option("header", "true").schema(schema_item).load(table_item_path)
    df_store_sales.registerTempTable("store_sales_table")
    df_item.registerTempTable("item_table")
    df_store_sales.cache().first()
    df_item.cache().first()

    val t0 = System.currentTimeMillis
    val lines = scala.io.Source.fromFile("/home/whassan/spark-sql-query-tests/src/main/scala/q26.sql").mkString
    val fin  = sqlContext.sql(lines)
    fin.collect()
    val t1 = System.currentTimeMillis
    println("****** Query 26 time(s) took: " + (t1 - t0).toFloat / 1000)
    fin.show()
    // println(fin.queryExecution.logical.numberedTreeString)
    // println("\n===================================\n")
    // println(fin.queryExecution.optimizedPlan.numberedTreeString)
    // println("\nExecuted Plan=====================\n")
    // println(fin.queryExecution.executedPlan.numberedTreeString)
    // println("\nSpark Plan=====================\n")
    // println(fin.queryExecution.sparkPlan.numberedTreeString)
    // println("\nStatistics=====================\n")
    // println(fin.queryExecution.analyzed.statistics.sizeInBytes)
    // println(df_item.queryExecution.analyzed.statistics.sizeInBytes)
    // println(df_store_sales.queryExecution.analyzed.statistics.sizeInBytes)
    println(fin.queryExecution.toString)
    println(":Done with Query 26")
  }
}
