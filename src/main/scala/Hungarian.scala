import breeze.linalg._

import scala.collection.immutable

  // Kuhn-Munkres typically stated with workers in rows and jobs in columns. Consider
  // transposing this for Breeze, since columns are dominant. Requires reorienting
  // workBoard and the order, (::, *) vs (::, *), of indexing matrix operations.

  // Copy cost matrix to mutable, markable cost slate:
  // - Transpose cost matrix as needed: shorter axis indexes worker; longer axis, jobs.
  // - Prepare sets of unmarked workers and jobs, so we can mark columns or rows containing 0.
  // - We subtract each worker's minimum job-cost from all the job costs for that worker.
  // - We subtract each job's minimum worker cost from all the worker costs for that job.
  // STEP 1. For each worker, SUBTRACT the min job cost from all the job costs for that worker.
  // For each worker, subtract their min job cost from all their job costs:
  // STEP 2. For each job, SUBTRACT the min worker cost from all the worker costs for that job.
  // (Note, for jobs with a worker cost set to 0 in Step 1, subtracting 0 obviously has no effect.)
  // For each job, subtract its min worker costs from all their worker costs:

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

class Hungarian(cost: DenseMatrix[Double]) {
  protected var costT: DenseMatrix[Double] =
    if (cost.rows <= cost.cols)
      cost.toDenseMatrix
    else
      cost.t.toDenseMatrix

  def reduceRowsCols(costT: DenseMatrix[Double]): DenseMatrix[Double] = {
    // Copy reoriented cost table
    val costX = costT.copy
    val workerMinJobCost = min(costX(::, *)).t
    costX(*, ::) :-= workerMinJobCost

    val jobMinWorkerCost = min(costX(*, ::))
    costX(::, *) :-= jobMinWorkerCost
    costX
  }

  case class Mark(worker: Boolean, index: Int)
  case class State(costX: DenseMatrix[Double],
                  // Default parameter values cover initialization
                   workerUnmarked: immutable.BitSet = immutable.BitSet(0 until costX.rows: _*),
                   jobUnmarked: immutable.BitSet = immutable.BitSet(0 until costX.cols: _*),
                   workerZeroJobs: DenseVector[immutable.BitSet] =
                     costX(*, ::).map(c => where(c :== 0.0)).
                       map(immutable.BitSet(_ : _*)),
                   jobZeroWorkers: DenseVector[immutable.BitSet] =
                     costX(::, *).map(c => where(c :== 0.0)).t.
                       map(immutable.BitSet(_ : _*))
                  ) {
    def getMark: Option[Mark] = {
      val bestWorkerMark: Int = workerUnmarked.maxBy(workerZeroJobs(_).size)
      val bestJobMark: Int = jobUnmarked.maxBy(jobZeroWorkers(_).size)
      val bestWorkerZeros = workerZeroJobs(bestWorkerMark).size
      val bestJobZeros = jobZeroWorkers(bestJobMark).size
      if (bestWorkerZeros == 0 && bestJobZeros == 0)
        None
      else if (bestJobZeros <= bestWorkerZeros)
        Some(Mark(worker = true, bestWorkerMark))
      else
        Some(Mark(worker = false, bestJobMark))
    }

    def markZeroIfAny(markMaybe : Option[Mark]): State = {
      markMaybe match {
        case None => this
        case Some(Mark(worker == true, bestWorkerMark)) => copy(
          costX = costX,
          workerUnmarked - bestWorkerMark,
          jobUnmarked,
          workerZeroJobs,
          jobZeroWorkers.map(_ - bestWorkerMark))
        case Some(Mark(worker == false, bestJobMark)) => copy(
          costX,
          workerUnmarked,
          jobUnmarked - bestJobMark,
          workerZeroJobs.map(_ - bestJobMark),
          jobZeroWorkers
        )
      }
    }
    def markAllZeros: State = {
      var state = copy()
      var markMaybe = getMark
      while (markMaybe.nonEmpty) {
        state = state.markZeroIfAny(markMaybe)
        markMaybe = state.getMark
      }
      state
    }

    def numWorkers: Int = costX.rows
    def numJobs: Int = costX.cols

    /** number of workers > Number of marks */
    def solved: Boolean = {
      numWorkers == (numWorkers - workerUnmarked.size) + (numJobs - jobUnmarked.size)
    }
    def reduceByMinUnmarked: State = {
      val minUnmarked = min(costX(workerUnmarked.toSeq, jobUnmarked.toSeq))
      // Subtract min unmarked cost from every unmarked worker;
      costX(workerUnmarked.toSeq, ::) :-= minUnmarked
      costX(::, 0 until numJobs filterNot(jobUnmarked(_))) :+= minUnmarked
      State(costX)  // Reset all marks.
    }
  }
  protected var costX: DenseMatrix[Double] = reduceRowsCols(costT)
  var state: State = State(costX)
  state = state.markAllZeros
  while (!state.solved) {
    state.reduceByMinUnmarked
    state = state.markAllZeros
  }

  // Fixme: Write code to read out solution, now that we have "minimum number of lines" = N.
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
