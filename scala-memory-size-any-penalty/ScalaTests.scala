import scala.language.postfixOps
//import SizeEstimator._

object ScalaTests {
  def int2(x : Array[Int]) = {
    for (i <- 0 to (x.length -1)) {
        x(i) = x(i) * 2
    }
  }

  def any2(x : Array[AnyVal]) = {
    for (i <- 0 to (x.length -1)) {
      x(i) = x(i).asInstanceOf[Int] * 2
    }
  }

  def test1() = {
    val tlen = 100000000
    var a = new Array[AnyVal](tlen)
    var b = new Array[Int](tlen)

    for (i <- 0 to (tlen - 1)) {
        a(i) = 1
        b(i) = 1
    }

    print("Size of a = ")
    val asize = SizeEstimator.estimate(a)
    print(asize)
    print(" per element = ")
    println(asize / tlen)

    print("Size of b = ")
    val bsize = SizeEstimator.estimate(b)
    print(bsize)
    print(" per element = ")
    println(bsize / tlen)

    val any_start = System.nanoTime()
    any2(a)
    val any_end   = System.nanoTime()
    val any_time = any_end.toDouble - any_start.toDouble

    val int_start = System.nanoTime()
    int2(b)
    val int_end   = System.nanoTime()
    val int_time = int_end.toDouble - int_start.toDouble

    val ret = any_time / int_time

    print("Multiply Any versus Int array by 2 multiplier = ")
    println(ret)
  }

  case class RowAny(x : Any, y : Any)
  case class RowDouble(x : Double, y : Double)

  def test2() = {
    val tlen = 10000000
    val ra : Array[RowAny] = for (i <- 1 to tlen toArray) yield RowAny(i.toDouble,i.toDouble)
    print("Size of ra = ")
    println(SizeEstimator.estimate(ra))

    val cstart = System.nanoTime()
    val rb : Array[RowDouble] = for (i <- ra) yield RowDouble(i.x.asInstanceOf[Double], i.y.asInstanceOf[Double])
    val cend   = System.nanoTime()
    val ctime = (cend.toDouble - cstart.toDouble) / rb.length.toDouble

    print("Size of rb = ")
    println(SizeEstimator.estimate(rb))

    print("Nano-second conversion time per element to Ret[Double, Double] = ")
    println(ctime)

    print("Total conversion time (s) to Ret[Double, Double] for ")
    print(tlen)
    print(" elements = ")
    println((cend.toDouble - cstart.toDouble) / 1000000000.0)
  }

  def main(args: Array[String]) {
    test1()
    test2()
  }
}
