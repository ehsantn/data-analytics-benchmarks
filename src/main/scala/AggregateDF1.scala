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
Example1: demonstrating a simple function that is not hard-coded
*/

object AggregateDF1 {
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
      .appName("AggregateDF1")
      .config("spark.sql.autoBroadcastJoinThreshold", "-1")
      .getOrCreate()

    import spark.implicits._
    val table_points_path = args(0)

    val schema_points = StructType(Array(
      StructField("id", LongType,true),
      StructField("x", DoubleType,true),
      StructField("y", DoubleType,true)))

    val df_points = spark.read.schema(schema_points).csv(table_points_path)

    df_points.registerTempTable("points")
    df_points.cache().first()
    // Starting time
    val t0 = System.currentTimeMillis
    val df1 = spark.sql("""SELECT
        id AS sid,
        sum(x)  AS sx,
        sum(y) AS sy
        FROM points 
        GROUP BY id
        """ )
    // head to intiate lazy evaluation
    df1.cache.head
    val t1 = System.currentTimeMillis
    // From spark website, there should be a good way
   // Measure time
    println("****** AggregateDF1 time(s) took: " + (t1 - t0).toFloat / 1000)
    println(":Done with AggregateDF1")
  }
}
