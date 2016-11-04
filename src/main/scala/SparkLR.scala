import java.util.Random

import scala.math.exp

import breeze.linalg._

import org.apache.spark.sql.SparkSession

/**
 * Logistic regression based classification.
 * Usage: SparkLR [slices]
 *
 */
object SparkLR {
  val N = 200000000  // Number of data points
  val D = 10   // Number of dimensions
  val R = 0.7  // Scaling factor
  val ITERATIONS = 10
  val rand = new Random(42)

  case class DataPoint(x: Vector[Double], y: Double)

  def generatePoint(i: Int): DataPoint = {
    val y = if (i % 2 == 0) -1 else 1
    val x = DenseVector.fill(D) {rand.nextGaussian + y * R}
    DataPoint(x, y)
  }

  def parseVector(line: String): DenseVector[Double] = {
    DenseVector(line.split(',').map(_.toDouble))
  }
  def main(args: Array[String]) {

    val spark = SparkSession
      .builder
      .appName("SparkLR")
      .getOrCreate()

    val sc = spark.sparkContext
    // val numSlices = if (args.length > 0) args(0).toInt else 2
    // val points = sc.parallelize(1 until N).map(generatePoint).cache()
    val lines = spark.read.textFile(System.getenv("HOME")+"/.julia/v0.5/HPAT/input_data/logistic_regression_64m.csv").rdd
    val points = lines.map(parseVector _).cache()
    points.first()

    val t0 = System.currentTimeMillis()
    // Initialize w to a random value
    var w = DenseVector.fill(D) {2 * rand.nextDouble - 1}
    println("Initial w: " + w)

    def grad(p: DenseVector[Double]) : DenseVector[Double] = {
        val x = p(0 to D-1)
        val y = p(D)
        val res = x * (1 / (1 + exp(-y * (w.dot(x)))) - 1) * y
        return res
    }
    for (i <- 1 to ITERATIONS) {
      println("On iteration " + i)
      val gradient = points.map(grad).reduce(_ + _)
      w -= gradient
    }

    val t1 = System.currentTimeMillis()
    println("Final w: " + w)
    val time = (t1-t0).toFloat/1000
    println("exec time "+(time)+"s")
    spark.stop()
  }
}
// scalastyle:on println
