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

object KernelScore3 {

  def main(args: Array[String]) {
    val spark = SparkSession
      .builder()
      .appName("KernelScore2")
      .getOrCreate()

    import spark.implicits._
    var n = 10e3.toInt
    if(args.length>0)
        n = args(0).toInt
    println("table size: ",n)

    var num_cores = 288
    if(args.length>1)
        num_cores = args(1).toInt
    println("total cores: ", num_cores)

    val rdd = spark.sparkContext.parallelize(1 until n, 2*num_cores).map(_ => Row(
        scala.util.Random.nextInt(256).toLong,Math.random(),Math.random()))
    val schema_points = StructType(Array(
        StructField("id", LongType,true),
        StructField("x", DoubleType,true),
        StructField("y", DoubleType,true)))
    val df_points = spark.createDataFrame(rdd, schema_points).repartition(2*num_cores)
    df_points.registerTempTable("points")
    df_points.cache().first()
    println("input read done")
    val b = 0.5 // bandwidth
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
        val minDis = dis.min
        var expSum = 0.0
        idx = 0
        while(idx<N) {
            expSum += math.exp(dis(idx)-minDis)
            idx += 1
        }
        minDis - math.log(b*N) + math.log(expSum)
    }
    
    spark.udf.register("scoreUDF", score _)
    // Starting time
    val t0 = System.currentTimeMillis
    val res = spark.sql("select sum(scoreUDF(x)) from points").collect()
    val t1 = System.currentTimeMillis
    // Measure time
    println("****** KernelScore3 time(s) took: " + (t1 - t0).toFloat / 1000)
    println(":Done with KernelScore3 ", res(0).get(0))
  }
}
