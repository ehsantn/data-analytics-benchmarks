build:	ScalaTests.class SizeEstimator.class

ScalaTests.class:	ScalaTests.scala SizeEstimator.class
	scalac -cp . ScalaTests.scala

SizeEstimator.class:	SizeEstimator.scala
	scalac SizeEstimator.scala

run:
	scala -J-Xms256m -J-Xmx4096m -cp . ScalaTests	
