build:	ScalaTests.class SizeEstimator.class

ScalaTests.class:	ScalaTests.scala SizeEstimator.class
	scalac -cp . ScalaTests.scala

SizeEstimator.class:	SizeEstimator.scala
	scalac SizeEstimator.scala

run:
	scala -J-Xms4096m -J-Xmx8096m -J-verbose:gc -J-XX:+PrintGCTimeStamps -J-XX:+PrintGCDetails -cp . ScalaTests
