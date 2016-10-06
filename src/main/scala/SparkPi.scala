// scalastyle:off println
import scala.math.random
import org.apache.spark.sql.SparkSession

/** Computes an approximation to pi */
object SparkPi {
  def main(args: Array[String]) {
    val spark = SparkSession
      .builder
      .appName("Spark Pi")
      .getOrCreate()
    //val slices = if (args.length > 0) args(0).toInt else 2
    val slices = 10000 
    val n = math.min(100000L * slices, Int.MaxValue).toInt // avoid overflow
    val t0 = System.currentTimeMillis
    val count = spark.sparkContext.parallelize(1 until n, slices).map { i =>
      val x = random * 2 - 1
      val y = random * 2 - 1
      if (x*x + y*y < 1) 1 else 0
    }.reduce(_ + _)
    println("Pi is roughly " + 4.0 * count / (n - 1))
    val t1 = System.currentTimeMillis
    println("** pi execution time(s) took: " + (t1 - t0).toFloat / 1000)
    spark.stop()
  }
}
// scalastyle:on println