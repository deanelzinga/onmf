package com.deanelzinga.onmf
class HungarianRoaring {

}

import org.roaringbitmap.RoaringBitmap
object HungarianRoaring {
  def main(args: Array[String]): Unit = {
    val roaring10 = RoaringBitmap.bitmapOf(0 to 10: _*)
    val roaring10even = RoaringBitmap.bitmapOf(0 to 10 by 2: _*)

    println(roaring10)
    println(roaring10even)
    println()
  }
}
