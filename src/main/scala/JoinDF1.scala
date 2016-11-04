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

object JoinDF1 {
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
      .appName("JoinDF1")
      .config("spark.sql.autoBroadcastJoinThreshold", "-1")
      .getOrCreate()

    import spark.implicits._
    val table_points_path = args(0)

    val schema_points1 = StructType(Array(
      StructField("id", LongType,true),
      StructField("x", DoubleType,true),
      StructField("y", DoubleType,true)))
    val schema_points2 = StructType(Array(
      StructField("id", LongType,true),
      StructField("x1", DoubleType,true),
      StructField("y1", DoubleType,true)))


    val df_points1 = spark.read.schema(schema_points1).csv(table_points_path)
    val df_points2 = spark.read.schema(schema_points2).csv(table_points_path)

    df_points1.registerTempTable("points1")
    df_points2.registerTempTable("points2")
    df_points1.cache().first()
    df_points2.cache().first()
    // Starting time
    val t0 = System.currentTimeMillis
    val df1 = spark.sql("""SELECT
        points1.id, x, y, x1, y1
        FROM points1 INNER JOIN points2 ON points1.id=points2.id
        """ )
    // head to intiate lazy evaluation
    df1.cache.head
    val t1 = System.currentTimeMillis
    // From spark website, there should be a good way
   // Measure time
    println("****** JoinDF1 time(s) took: " + (t1 - t0).toFloat / 1000)
    println(":Done with JoinDF1")
  }
}
