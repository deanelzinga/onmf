package com.deanelzinga.gopa

import breeze.linalg.DenseMatrix
import org.scalatest.FunSuite
//
// import kuhnmunkres.HungarianT
class KuhnMunkresTest extends FunSuite {

  // Tests go here
  test("State(DenseMatrix.ones[Double](1,1)).toString should be \"       0\n  0    1\n\"") {

    assert(new HungarianT(DenseMatrix.ones[Double](1, 1)).
      ZeroMarking(DenseMatrix.ones[Double](1, 1)).toString ==
      "       0\n  0    1\n")
  }

}

