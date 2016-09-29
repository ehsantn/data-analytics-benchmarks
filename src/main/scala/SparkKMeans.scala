
// scalastyle:off println

import breeze.linalg.{squaredDistance, DenseVector, Vector}

import org.apache.spark.sql.SparkSession

/**
 * K-means clustering.
 *
 */
object SparkKMeans {

  def parseVector(line: String): Vector[Double] = {
    DenseVector(line.split(',').map(_.toDouble))
  }

  def closestPoint(p: Vector[Double], centers: Array[Vector[Double]]): Int = {
    var bestIndex = 0
    var closest = Double.PositiveInfinity

    for (i <- 0 until centers.length) {
      val tempDist = squaredDistance(p, centers(i))
      if (tempDist < closest) {
        closest = tempDist
        bestIndex = i
      }
    }

    bestIndex
  }

  def main(args: Array[String]) {

//    if (args.length < 3) {
  //    System.err.println("Usage: SparkKMeans <file> <k> <convergeDist>")
//      System.exit(1)
//    }

    val spark = SparkSession
      .builder
      .appName("SparkKMeans")
      .getOrCreate()

    //val lines = spark.read.textFile(args(0)).rdd
    //val data = lines.map(parseVector _).cache()
    //val K = args(1).toInt
    //val convergeDist = args(2).toDouble
    val lines = spark.read.textFile(System.getenv("HOME")+"/.julia/v0.5/HPAT/input_data/kmeans_large.csv").rdd
    val data = lines.map(parseVector _).cache()
    val K = 5
    val iterations = 10

    val t0 = System.currentTimeMillis
    
    val kPoints = data.takeSample(withReplacement = false, K, 1)

    for( i <- 1 to iterations) {
      val closest = data.map (p => (closestPoint(p, kPoints), (p, 1)))

      val pointStats = closest.reduceByKey{case ((p1, c1), (p2, c2)) => (p1 + p2, c1 + c2)}

      val newPoints = pointStats.map {pair =>
        (pair._1, pair._2._1 * (1.0 / pair._2._2))}.collectAsMap()

      for (newP <- newPoints) {
        kPoints(newP._1) = newP._2
      }
      println("Finished iteration " + i )
    }

    println("Final centers:")
    kPoints.foreach(println)
    val t1 = System.currentTimeMillis
    println("*** Kmeans execution time(s) took: " + (t1 - t0).toFloat / 1000)
    spark.stop()
  }
}
// scalastyle:on println
