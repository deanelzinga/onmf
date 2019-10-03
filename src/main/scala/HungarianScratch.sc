import breeze.linalg._
import breeze.numerics._
import scala.collection.mutable
import org.roaringbitmap.RoaringBitmap


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

cost2r = cost2r(::, *) -  min(cost2r(*, ::))
argmin(cost2r)
mutable.BitSet((0 until 10) : _*)
RoaringBitmap.bitmapOf((0 until 10) : _*)

val vec25 = DenseVector.range(0, 25)
val indices = mutable.BitSet.empty ++ (0 until 10)

vec25.mapPairs((k,v) => (if (indices(k)) v else None))
var a = 2
var b = 3
var(c, d) = (3, 2)

1 to 3 flatMap(i => 5 to 8 map(j => (i, j)))

(1 to 10 toSeq) -- (1 to 3 toSet)