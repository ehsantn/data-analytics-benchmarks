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
import scala.language.reflectiveCalls
import java.lang.management.ManagementFactory
import scala.collection.JavaConversions._
/**
*/

object KernelScore2 {
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
      .appName("KernelScore2")
      .config("spark.sql.autoBroadcastJoinThreshold", "-1")
      .getOrCreate()

    import spark.implicits._
    val table_points_path = args(0)
    val b = 3.0 // bandwidth
    val points = Array(-1.0, 2.0, 5.0)
    val N = points.length
    def score(x: Double): Double = {
        val dis = new Array[Double](N)
        var idx = 0
        while(idx<N) {
            val z_k = points(idx)
            dis(idx) = -(x - z_k) * (x - z_k) / (2*b*b)
            idx += 1
        }
        //val dis = points.map { z_k =>
        //    -(x - z_k) * (x - z_k) / (2*b*b)
        // }
        val minDis = dis.min
        var expSum = 0.0
        idx = 0
        while(idx<N) {
            expSum += math.exp(dis(idx)-minDis)
            idx += 1
        }
        //val exps = dis.map(d => math.exp(d-minDis))
        //minDis - math.log(b*N) + math.log(exps.sum)
        minDis - math.log(b*N) + math.log(expSum)
    }
    val schema_points = StructType(Array(
        StructField("id", LongType,true),
        StructField("x", DoubleType,true),
        StructField("y", DoubleType,true)))

    val df_points = spark.read.schema(schema_points).csv(table_points_path)
    df_points.registerTempTable("points")
    df_points.cache().first()
    spark.udf.register("scoreUDF", score _)
    // Starting time
    val t0 = System.currentTimeMillis
    val res = spark.sql("select sum(scoreUDF(x)) from points").collect()
    val t1 = System.currentTimeMillis
    // From spark website, there should be a good way
    // Measure time
    println("****** KernelScore2 time(s) took: " + (t1 - t0).toFloat / 1000)
    println(":Done with KernelScore2 ", res)
  }
}
