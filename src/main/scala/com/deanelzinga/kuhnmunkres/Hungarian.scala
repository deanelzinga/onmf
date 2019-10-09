package com.deanelzinga.kuhnmunkres

import breeze.linalg._
import com.deanelzinga.kuhnmunkres.Hungarian._

import scala.collection.mutable

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
  protected var costT: DenseMatrix[Double] = transposeIfAny(cost)
  protected var costX: DenseMatrix[Double] = costT.copy
  def numWorkers: Int = costX.rows
  def numJobs: Int = costX.cols
  reduceRows(costX)
  reduceCols(costX)
  private def newWorkerUnmarked() = mutable.BitSet(0 until numWorkers: _*)
  private def newJobUnmarked() = mutable.BitSet(0 until numJobs: _*)
  private def newWorkerZeroJobs(): DenseVector[mutable.BitSet] = costX(*, ::).
    map(c => where(c :== 0.0)).
    map(mutable.BitSet(_: _*))
  private def newJobZeroWorkers(): DenseVector[mutable.BitSet] = costX(::, *).
    map(c => where(c :== 0.0)).t.
    map(mutable.BitSet(_: _*))
  case class Mark(worker: Boolean, index: Int)

  case class State(costX: DenseMatrix[Double] = costX,
                   // Default parameter values cover initialization
                   var workerUnmarked: mutable.BitSet = newWorkerUnmarked(),
                   var jobUnmarked: mutable.BitSet = newJobUnmarked(),
                   var workerZeroJobs: DenseVector[mutable.BitSet] = newWorkerZeroJobs(),
                   var jobZeroWorkers: DenseVector[mutable.BitSet] = newJobZeroWorkers()) {
    //def numWorkers = costX.rows
    //protected def numjobs = costX.cols
    // def workerNumZeros(worker: Int): Int = workerZeroJobs(worker).size
    def reset(): Unit = {
      workerUnmarked = newWorkerUnmarked()
      jobUnmarked = newJobUnmarked()
      workerZeroJobs = newWorkerZeroJobs()
      jobZeroWorkers = newJobZeroWorkers()
    }

    override def toString: String = {
      val jobIndices = (0 until numJobs) map( job => f"$job%4d") mkString("    "+" job", "", "\n")
      val jobZeroSizes = (0 until numJobs) map(
        job => {
          val numZeroWorkers = jobZeroWorkers(job).size
          if (jobUnmarked(job)) "   |" else f"$numZeroWorkers%4d"
        }) mkString(" wrk"+" mrk", "", "\n")
      val lines = (0 until numWorkers) map {
        worker => {
          val numZeroJobs = workerZeroJobs(worker).size
          Seq(f"$worker%4d",
            if (workerUnmarked(worker)) "  --" else f"$numZeroJobs%4d",
            costX(worker, ::).t.toScalaVector().map(c => if (c == 0.0) "   ." else f"$c%4.0f").mkString
          ).mkString
        }
      }
      Seq(jobIndices, jobZeroSizes, lines.mkString("","\n","\n")).mkString
    }

    def getMarkIfAny: Option[Mark] = {
      var bestWorkerMark = -1
      var bestWorkerZeros = 0
      var bestJobMark = -1
      var bestJobZeros = 0
      if (workerUnmarked.nonEmpty) {
        bestWorkerMark = workerUnmarked.maxBy[Int](workerZeroJobs(_).size)
        bestWorkerZeros = workerZeroJobs(bestWorkerMark).size
      }
      if (jobUnmarked.nonEmpty) {
        bestJobMark = jobUnmarked.maxBy[Int](jobZeroWorkers(_).size)
        bestJobZeros = jobZeroWorkers(bestJobMark).size
      }
      if (bestWorkerZeros == 0 && bestJobZeros == 0)
        None
      else if (bestJobZeros <= bestWorkerZeros)
        Some(Mark(worker = true, bestWorkerMark))
      else
        Some(Mark(worker = false, bestJobMark))
    }

    def markZeroIfAny(markMaybe : Option[Mark]): Unit = {
      markMaybe match {
        case None => ()
        case Some(Mark(true, bestWorkerMark)) =>
          workerUnmarked -= bestWorkerMark
          workerZeroJobs(bestWorkerMark).foreach(job => jobZeroWorkers(job) -= bestWorkerMark)
        case Some(Mark(false, bestJobMark)) =>
          jobUnmarked -= bestJobMark
          jobZeroWorkers(bestJobMark).foreach(worker => workerZeroJobs(worker) -= bestJobMark)
        case _ => throw new RuntimeException("Unhandled markMaybe case")
      }
    }

    def markAllZeros(): Unit = {
      var markMaybe = getMarkIfAny
      while (markMaybe.nonEmpty) {
        state.markZeroIfAny(markMaybe)
        markMaybe = state.getMarkIfAny
      }
    }

    /** number of workers > Number of marks */
    def solved: Boolean = {
      numWorkers == (numWorkers - workerUnmarked.size) + (numJobs - jobUnmarked.size)
    }

    def reduceByMinUnmarked(): Unit = {
      val minUnmarked =
        // min(costX(workerUnmarked.toSeq, jobUnmarked.toSeq))
        min(costX(new BitVector(data = workerUnmarked),
          new BitVector(data = jobUnmarked)))
      // Subtract min unmarked cost from every unmarked worker;
      costX(workerUnmarked.toSeq, ::) :-= minUnmarked
      // Add min unmarked cost to every MARKED column:
      costX(::, 0 until numJobs filterNot(jobUnmarked(_))) :+= minUnmarked
      reset() // Reset all marks.
    }
  }

  object State {

  }

  def solve(): Unit = {
    while (!state.solved) {
      state.reduceByMinUnmarked()
      state.markAllZeros()
    }
  }

  var state: State = State()
  state.markAllZeros()

  // Fixme: Write code to read out solution, now that we have "minimum number of lines" = N.
}

object Hungarian {

  def apply(cost: DenseMatrix[Double]): Hungarian = {
    new Hungarian(cost)
  }
  def transposeIfAny(cost: DenseMatrix[Double]): DenseMatrix[Double] = {
    if (cost.rows <= cost.cols)
      cost.toDenseMatrix
    else
      cost.t.toDenseMatrix
  }
  def reduceRows(costX: DenseMatrix[Double]): Unit = {
    val jobMinWorkerCost = min(costX(*, ::))
    costX(::, *) :-= jobMinWorkerCost
  }

  def reduceCols(costX: DenseMatrix[Double]): Unit = {
    val workerMinJobCost = min(costX(::, *)).t
    costX(*, ::) :-= workerMinJobCost
  }

  /**
   *
   */
  def main(args: Array[String]): Unit = {
    val m4by4 = DenseMatrix((1 to 16).map(_.toDouble).toArray).reshape(4,4).t
    val m4by4sqr = m4by4 *:* m4by4
    val h = Hungarian(cost = m4by4sqr)
    println(h.state.toString)
    h.state.reduceByMinUnmarked()
  }
}
