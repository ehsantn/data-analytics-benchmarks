import scala.language.existentials
/**
*/

object KernelScoreS {

  def main(args: Array[String]) {
    var n = 10e3.toInt
    if(args.length>0)
        n = args(0).toInt
    println("table size: ",n)

    var X = new Array[Double](n)
    for (i <- 0 to n-1){
        X(i) = Math.random();
    }
    println("input generate done")
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
    
    var allsum = 0.0
    // Starting time
    val t0 = System.currentTimeMillis
    for (i <- 0 to n-1){
        allsum += score(X(i))
    }
    val t1 = System.currentTimeMillis
    // Measure time
    println("****** KernelScoreS time(s) took: " + (t1 - t0).toFloat / 1000)
    println(":Done with KernelScoreS ", allsum)
  }
}
