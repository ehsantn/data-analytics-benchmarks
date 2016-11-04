import scala.math.random
import org.apache.spark.sql.SparkSession

object OneDSumFilter {
  
  def parseVector(line: String): Double = {
    line.toDouble
  }

  def main(args: Array[String]) {
    val spark = SparkSession
      .builder
      .appName("Spark Pi")
      .getOrCreate()
    val lines = spark.read.textFile(System.getenv("HOME")+"/.julia/v0.5/HPAT/input_data/1D_large.csv").rdd
    val data = lines.map(parseVector _).cache()
    data.first()
    
    val t0 = System.currentTimeMillis
    val res = data.map(p=> if (p>0.8) p else 0.0).reduce(_ + _)
    println("sum is  " + res)
    val t1 = System.currentTimeMillis
    println("** sum filter execution time(s) took: " + (t1 - t0).toFloat / 1000)
    spark.stop()
  }
}
