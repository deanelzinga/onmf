Global / onChangedBuildSource := ReloadOnSourceChanges
name := "kuhnmunkres"
version := "0.1"
scalaVersion := "2.13.1"
libraryDependencies += "org.scalactic" %% "scalactic" % "3.0.8"
libraryDependencies += "org.scalatest" %% "scalatest" % "3.0.8" % "test"
// resolvers += "Artima Maven Repository" at "https://repo.artima.com/releases"
//logBuffered in Test := false

libraryDependencies  ++= Seq(
  // Last stable release
  "org.scalanlp" %% "breeze" % "1.0",

  // Native libraries are not included by default. add this if you want them
  // Native libraries greatly improve performance, but increase jar sizes.
  // It also packages various blas implementations, which have licenses that may or may not
  // be compatible with the Apache License. No GPL code, as best I know.
  "org.scalanlp" %% "breeze-natives" % "1.0",

  // The visualization library is distributed separately as well.
  // It depends on LGPL code
  "org.scalanlp" %% "breeze-viz" % "1.0"
)
libraryDependencies += "org.roaringbitmap" % "RoaringBitmap" % "0.8.11"
libraryDependencies += "org.openjdk.jol" % "jol-core" % "0.9"

// Logging, from https://github.com/lightbend/scala-logging
libraryDependencies += "com.typesafe.scala-logging" %% "scala-logging" % "3.9.2"
libraryDependencies += "ch.qos.logback" % "logback-classic" % "1.2.3"
libraryDependencies += "it.unimi.dsi" % "fastutil" % "8.3.0"