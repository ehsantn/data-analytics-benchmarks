name := "Query26"
version := "0.1"
scalaVersion := "2.11.8"
libraryDependencies += "org.apache.spark" %% "spark-core" % "2.0.0"
libraryDependencies += "org.apache.spark" %% "spark-sql" % "2.0.0"
libraryDependencies += "org.apache.spark" %% "spark-mllib" % "2.0.0"

libraryDependencies += "org.scalanlp" %% "breeze" % "0.11.2"
//libraryDependencies  ++= Seq(
    // other dependencies here
//    "org.scalanlp" %% "breeze" % "0.12",
    // native libraries are not included by default. add this if you want them (as of 0.7)
    // native libraries greatly improve performance, but increase jar sizes. 
    // It also packages various blas implementations, which have licenses that may or may not
    // be compatible with the Apache License. No GPL code, as best I know.
  //  "org.scalanlp" %% "breeze-natives" % "0.12",
    // the visualization library is distributed separately as well. 
    // It depends on LGPL code.
  //  "org.scalanlp" %% "breeze-viz" % "0.12"
//)
retrieveManaged := true
