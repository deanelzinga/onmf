import breeze.linalg._
import breeze.numerics._
import scala.collection.mutable
import org.roaringbitmap.RoaringBitmap
import com.deanelzinga.gopa._
import scala.collection.immutable

val mtri0 = DenseMatrix((0.0, 0.0, 1.0, 1.0, 1.0),
  (0.0, 0.0, 0.0, 1.0, 1.0),
  (1.0, 0.0, 0.0, 0.0, 1.0),
  (1.0, 1.0, 0.0, 0.0, 0.0),
  (1.0, 1.0, 1.0, 0.0, 0.0))
val cost = DenseMatrix.rand[Double](4, 4)
val working = cost.copy

// For each worker, subtract their min job cost from all their job costs:
val minJobCostPerWorker = min(working(*, ::))
working(::, *) :-= minJobCostPerWorker

// For each job, subtract its min worker costs from all their worker costs:
val minWorkerCostPerJob = min(working(::, *)).t
working(*, ::) :-= minWorkerCostPerJob

// This leaves zeros in the rows and columns

val jobCostZeros = working(::, *).
  map(jobCost => where(jobCost :== 0.0)).
  t.
  map(_.toSet[Int])
val workerCostZeros = working(*, ::).
  map(workerCost => where(workerCost :== 0.0)).
  map(_.toSet[Int])
val jobCostMax = max(jobCostZeros.map(_.size))
val workerCostMax = max(workerCostZeros.map(_.size))

val cost2 = DenseMatrix.rand[Double](5, 4)
var cost2r = cost2.t.toDenseMatrix
min(cost2r(::, *))
cost2r = cost2r(*, ::) - min(cost2r(::, *)).t

cost2r = cost2r(::, *) - min(cost2r(*, ::))
argmin(cost2r)
mutable.BitSet((0 until 10) : _*)
RoaringBitmap.bitmapOf((0 until 10) : _*)

val vec25 = DenseVector.range(0, 25)
val indices = mutable.BitSet.empty ++ (0 until 10)

vec25.mapPairs((k,v) => if (indices(k)) v else None)
var a = 2
var b = 3
var(c, d) = (3, 2)

(1 to 3).flatMap(i => 5 to 8 map(j => (i, j)))

Matrix.ones[Double](3,3)
val m0 = DenseMatrix(1 to 16 map(_.toDouble)).reshape(4,4).t
val m1 = m0 *:* m0
val m = DenseMatrix((1.0, 2.0, 3.0),(4.0, 5.0, 6.0), (7.0, 8.0, 9.0))
//val hs = "----\n" + h.state.toString + "----"
val m0112 = DenseMatrix((0.0, 1.0), (1.0, 2.0))
val mScratch = m1.copy
HungarianT.reduceCols(mScratch)
println(mScratch)
HungarianT.reduceRows(mScratch)
println(mScratch)
val h = new HungarianT(m1)
h.zeroMarking.toString
h.zeroMarking.relaxByMinUnmarked()
h.zeroMarking.toString
println(h.zeroMarking.toString)
h.zeroMarking.minMarkZeros()
h.zeroMarking.relaxByMinUnmarked()
val d1to16 = DenseVector((1 to 16).map(_.toDouble).toArray)

val bits: BitVector = DenseVector((1 to 16).map(_.toDouble).toArray) <:< 8.0

d1to16(bits)
