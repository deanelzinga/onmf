import breeze.linalg._

import scala.collection.{SortedSet, mutable}

//import scala.collection.mutable
import scala.collection.immutable
import org.roaringbitmap.RoaringBitmap
import scala.collection.mutable.BitSet

class Hungarian(cost: DenseMatrix[Double]) {
  // Kuhn-Munkres typically stated with workers in rows and jobs in columns. Consider
  // transposing this for Breeze, since columns are dominant. Requires reorienting
  // workBoard and the order, (::, *) vs (::, *), of indexing matrix operations.

  // Copy cost matrix to mutable, markable cost slate:
  // - Transpose cost matrix as needed: shorter axis indexes worker; longer axis, jobs.
  // - Prepare sets of unmarked workers and jobs, so we can mark columns or rows containing 0.
  // - We subtract each worker's minimum job-cost from all the job costs for that worker.
  // - We subtract each job's minimum worker cost from all the worker costs for that job.

  protected var costx: DenseMatrix[Double] =
    if (cost.rows <= cost.cols)
      cost.toDenseMatrix
    else
      cost.t.toDenseMatrix
  protected var numWorkers: Int = costx.rows
  protected var numJobs: Int = costx.cols

  // STEP 1. For each worker, subtract the min job cost from all the job costs for that worker.
  // For each worker, subtract their min job cost from all their job costs:
  private val minJobCostPerWorker = min(costx(::, *)).t
  costx(*, ::) :-= minJobCostPerWorker

  // STEP 2. For each job, subtract the min worker cost from all the worker costs for that job.
  // (Note, for jobs with a worker cost set to 0 in Step 1, subtracting 0 obviously has no effect.)
  // For each job, subtract its min worker costs from all their worker costs:
  private val minWorkerCostPerJob = min(costx(*, ::))
  costx(::, *) :-= minWorkerCostPerJob

  // STEP 3. Starting from an unmarked, transformed cost matrix...
  // - Mark all the 0s of the transformed cost matrix with a minimum number of marks on
  // workers or jobs (rows or columns).

  // Find and make next mark, if any:
  // At each step, we use a greedy algorithm:
  // - Find the max count of unmarked, 0-cost jobs for any unmarked worker.
  // - Find the maximum number of unmarked, 0-cost workers for any unmarked job.
  // - If there are no more unmarked, 0-cost jobs, go to STEP 4.
  // - If either is greater, choose it. If they are equal, prefer the worker.
  // - Mark that worker or job (remove it from unmarked workers or jobs).
  // - Remove that worker or job from all jobs' or workers' sets of 0-cost workers or jobs.

  // STEP 4. TEST FOR OPTIMALITY: (i) If the minimum of 0-covering marks is N,
  // the number of workers, then an optimal assignment of 0s is possible and "we are
  // finished"-- Go to Step 6.

  // STEP 5. Find the minimum cost entry still unmarked.
  // - Subtract this entry from each unmarked worker.
  // - Add this entry to each MARKED job. Return to STEP 3.

  // STEP 6: Read off an available assignment from the covered 0s (not completely trivial).
  var workerZeroJobs: DenseVector[immutable.BitSet] = DenseVector(immutable.BitSet.empty)

  var jobZeroWorkers: DenseVector[immutable.BitSet] = DenseVector(immutable.BitSet.empty)
  var jobUnmarked: immutable.BitSet = immutable.BitSet.empty
  var workerUnmarked: immutable.BitSet = immutable.BitSet.empty
  case class Mark(worker: Boolean, index: Int)
  def getMark(workerUnmarked: immutable.BitSet,
              workerZeroJobs: DenseVector[immutable.BitSet],
              jobUnmarked: immutable.BitSet,
              jobZeroWorkers: DenseVector[immutable.BitSet]
               ): Option[Mark] = {
    var bestWorkerMark: Int = workerUnmarked.maxBy(workerZeroJobs(_).size)
    var bestJobMark: Int = jobUnmarked.maxBy(jobZeroWorkers(_).size)
    var bestWorkerZeros = workerZeroJobs(bestWorkerMark).size
    var bestJobZeros = jobZeroWorkers(bestJobMark).size
    if (bestWorkerZeros == 0 && bestJobZeros == 0)
      None
    else if (bestJobZeros <= bestWorkerZeros)
      Some(Mark(worker = true, bestWorkerMark))
    else
      Some(Mark(worker = false, bestJobMark))
  }
  // For each worker (::), over all jobs (*), map out set of indices where cost is zero:
  workerZeroJobs = costx(::, *).map(jobCost => where(jobCost :== 0.0)).
    t.map(immutable.BitSet(_: _*))
  // For each job (::), over all workers (*),
  jobZeroWorkers = costx(*, ::).map(workerCost => where(workerCost :== 0.0)).
    map(immutable.BitSet(_: _*))
  jobUnmarked = immutable.BitSet((0 until numJobs): _*)
  workerUnmarked = immutable.BitSet((0 until numWorkers): _*)
  var mark: Option[Mark] = getMark(workerUnmarked, workerZeroJobs, jobUnmarked, jobZeroWorkers)
  while (mark.nonEmpty) {
    if (mark.get.worker) {
      workerUnmarked = workerUnmarked - mark.get.index
      jobZeroWorkers = jobZeroWorkers.map(_ - mark.get.index)
    } else {
      jobUnmarked = jobUnmarked - mark.get.index
      workerZeroJobs = workerZeroJobs.map(_ - mark.get.index)
    }
    mark = getMark(workerUnmarked, workerZeroJobs, jobUnmarked, jobZeroWorkers)
  }
  while ((numWorkers - workerUnmarked.size) + (numJobs - jobUnmarked.size) < numWorkers) {
    val cross = workerUnmarked.flatMap(w => jobUnmarked.map(j => (w, j)))
    val wjMin = cross.minBy(wj => costx(wj._1, wj._2))
    val minUnmarkedCost = costx(wjMin._1, wjMin._2)
    // Subtract min unmarked cost from every unmarked worker;
    costx(workerUnmarked.toSeq, ::) :-= minUnmarkedCost
    costx(::, 0 until numJobs filterNot(jobUnmarked(_))) :+= minUnmarkedCost
    // FIXME: Add steps: clear marks, redo marks.
  }
}

object Hungarian {
  /**
   *
   * @return Vector of column indices for lowest cost.
   */

  def optimize(cost: Matrix[Double]): Vector[Int] = {
    val rows = cost.rows
    val c = cost.copy
    // For every row, subtract the minimum of that row

    val rowRange = DenseVector.range(0, rows)
    rowRange
  }
}
