package com.deanelzinga.onmf

import breeze.linalg._
import it.unimi.dsi.fastutil.ints.{IntComparator, IntHeapIndirectPriorityQueue}

import scala.collection.mutable

/**
 * Kuhn-Munkres typically stated with workers in rows and jobs in columns. This is transposed, because
 * Breeze (and BLAS) operates on column-major matrices. In other words:
 * - Columns are the primary objects in a Breeze (and BLAS) matrix;
 * - Workers are the primary objects in Kuhn-Munkres;
 * - It makes some sense, especially in cases where numWorkers << numJobs, for each worker to be a vector
 * along the longer axis, the axis of jobs.
 * Example: Cost matrix is 10 x 1000:
 * - Each worker is represented along the short axis (10 workers).
 * - Each worker has costs associated with each job (1000 jobs).
 * - For a column-major system like Breeze (and BLAS), it makes sense to represent this as 1000 rows (jobs) by
 * 10 column-vectors (workers).
 *
 * For all that, this is currently unused and not directly needed, built on speculation it might be needed for Gopa.
 */
@deprecated
class HungarianT(cost: DenseMatrix[Double]) {
  // Originally the algorithm was stated with rows fewer than columns. This matters for some of the order
  // of operations during the algorithm. For Breeze, which is column oriented, it probably makes more sense
  // for columns to be longer than rows. In anticipation of this change, the shorter axis is called workers
  // and the longer axis is called jobs, to separate the algorithmic decisionmaking from the implementation
  // details.
  protected var costT: DenseMatrix[Double] = {
    if (cost.cols <= cost.rows)
      cost.toDenseMatrix
    else
      cost.t.toDenseMatrix
  }

  protected var costX: DenseMatrix[Double] = costT.copy
  def reduceCols(): Unit = {  // Reduce columns
    val workerMinJobCost = min(costX(::, *)).t
    costX(*, ::) :-= workerMinJobCost
  }
  reduceCols()

  def reduceRows(): Unit = { // Reduce rows
    val jobMinWorkerCost = min(costX(*, ::))
    costX(::, *) :-= jobMinWorkerCost
  }
  reduceRows()
  def numWorkers: Int = costX.cols
  def numJobs: Int = costX.rows
  private def newWorkerUnmarked() = mutable.BitSet(0 until numWorkers: _*)
  private def newJobUnmarked() = mutable.BitSet(0 until numJobs: _*)

  private def newWorkerZeroJobs(): DenseVector[mutable.BitSet] = { costX(::, *).
    map(c => where(c :== 0.0)).t.
    map(mutable.BitSet(_: _*))
  }
  private def newJobZeroWorkers(): DenseVector[mutable.BitSet] =  { costX(*, ::).
    map(c => where(c :== 0.0)).
    map(mutable.BitSet(_: _*))
  }
  case class Mark(worker: Boolean, index: Int)

  /**
   *
   * @param costX transformed cost matrix, with workers the shortest axis
   * @param workerUnmarked
   * @param jobUnmarked
   * @param workerZeroJobs
   * @param jobZeroWorkers
   */
  case class ZeroMarking(costX: DenseMatrix[Double] = costX,
                         // Default parameter values cover initialization
                         var workerUnmarked: mutable.BitSet = newWorkerUnmarked(),
                         var jobUnmarked: mutable.BitSet = newJobUnmarked(),
                         var workerZeroJobs: DenseVector[mutable.BitSet] = newWorkerZeroJobs(),
                         var jobZeroWorkers: DenseVector[mutable.BitSet] = newJobZeroWorkers()) {

    /**
     * Resets the state after a reduction step:
     * - Resets all row and column marks along both axes.
     * - Recomputes sets of zero indices at every index; for every row,
     * the set of zero columns; for every column, the set of zero rows.
     */
    def resetMarks(): Unit = {
      workerUnmarked = newWorkerUnmarked()
      jobUnmarked = newJobUnmarked()
      workerZeroJobs = newWorkerZeroJobs()
      jobZeroWorkers = newJobZeroWorkers()
    }

    /**
     *
     * @return
     */
    override def toString: String = {
      val workerIndices = (0 until numWorkers) map (worker => f"$worker%4d") mkString("    " + " wrk", "", "\n")
      val workerZeroSizes = (0 until numWorkers) map (
        worker => {
          val numZeroJobs = workerZeroJobs(worker).size
          if (workerUnmarked(worker)) "   |" else f"$numZeroJobs%4d"
        }) mkString(" job" + " mrk", "", "\n")
      val lines = (0 until numJobs) map {
        job => {
          val numZeroWorkers = jobZeroWorkers(job).size
          Seq(f"$job%4d",
            if (jobUnmarked(job)) "  --" else f"$numZeroWorkers%4d",
            costX(job, ::).t.toScalaVector().map(c => if (c == 0.0) "   ." else f"$c%4.0f").mkString
          ).mkString
        }
      }
      Seq(workerIndices, workerZeroSizes, lines.mkString("", "\n", "\n")).mkString
    }

    /**
     * For either axis, 3 return values from a single O(N) search for the best next mark.
     *
     * @param axisMaxMark
     * @param axisNumUnmarkedZeros
     * @param axisMaxNumZeros
     */
    case class AxisMark(axisMaxMark: Int, axisMaxNumZeros: Int, axisNumUnmarkedZeros: Int)

    /**
     * Low-level, multi-return max function along either worker or job axis.
     *
     * @param axisUnmarked
     * @param axisZeroIndices
     * @return Optional 3 values in AxisMark: index with max number of zero indices available to mark,
     *         the number of zeros to mark at that index, and the number of unmarked indices with any zeros.
     */
    def nextAxisMarkOption(axisUnmarked: mutable.BitSet,
                           axisZeroIndices: DenseVector[mutable.BitSet]): Option[AxisMark] = {
      var axisMaxMark = -1
      var axisMaxNumZeros = 0
      var axisNumUnmarkedZeros = 0
      if (axisUnmarked.nonEmpty) {
        for (i <- axisUnmarked) {
          val iZeros = axisZeroIndices(i).size
          if (iZeros > 0) {
            axisNumUnmarkedZeros += 1
            if (iZeros > axisMaxNumZeros) {
              axisMaxMark = i
              axisMaxNumZeros = iZeros
            }
          }
        }
        if (axisNumUnmarkedZeros > 0)
          Some(AxisMark(axisMaxMark = axisMaxMark,
            axisNumUnmarkedZeros = axisNumUnmarkedZeros,
            axisMaxNumZeros = axisMaxNumZeros))
        else None
      }
      else None
    }

    /**
     *
     * @return
     */
    def nextMarkOption: Option[Mark] = {
      val workerMarkOption = nextAxisMarkOption(workerUnmarked, workerZeroJobs)
      val jobMarkOption = nextAxisMarkOption(jobUnmarked, jobZeroWorkers)
      (workerMarkOption, jobMarkOption) match {
        case (None, None) => None
        case (_, None) => throw new IllegalStateException("Unmarked zeros on job axis, not on worker axis.")
        case (None, _) => throw new IllegalStateException("Unmarked zeros on worker axis, not on job axis.")
        case (Some(AxisMark(workerMaxMark, workerMaxNumZeros, workerNumUnmarkedZeros)),
        Some(AxisMark(jobMaxMark, jobMaxNumZeros, jobNumUnmarkedZeros))) =>

          // If we can mark strictly more zeros on either axis, do it:
          if (jobMaxNumZeros < workerMaxNumZeros)
            Some(Mark(worker = true, workerMaxMark))
          else if (jobMaxNumZeros > workerMaxNumZeros)
            Some(Mark(worker = false, jobMaxMark))

          // Else choose shorter queue of unmarked rows or cols with zeros; for ties, choose worker (shorter axis).
          else if (jobNumUnmarkedZeros >= workerNumUnmarkedZeros)
            Some(Mark(worker = true, workerMaxMark))
          else
            Some(Mark(worker = false, jobMaxMark))
      }
    }

    def markZeroIfAny(markMaybe: Option[Mark]): Unit = {
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

    def minMarkZeros(): Unit = {
      var markMaybe = nextMarkOption
      while (markMaybe.nonEmpty) {
        zeroMarking.markZeroIfAny(markMaybe)
        markMaybe = zeroMarking.nextMarkOption
      }
    }

    /** number of workers > Number of marks */
    def solved: Boolean = {
      numWorkers == (numWorkers - workerUnmarked.size) + (numJobs - jobUnmarked.size)
    }

    def relaxByMinUnmarked(): Unit = {
      if (workerUnmarked.nonEmpty && jobUnmarked.nonEmpty) {
        val minUnmarked = min(costX(jobUnmarked.toSeq, workerUnmarked.toSeq))

        // Subtract min unmarked cost from every unmarked worker;
        costX(::, workerUnmarked.toSeq) :-= minUnmarked

        // Add min unmarked cost to every MARKED column:
        costX((0 until numJobs).filterNot(jobUnmarked(_)), ::) :+= minUnmarked
      }
      else ()
    }

    // FIXME: Refactor HungarianT.solve() to be HungarianT.state.solve()
    def relaxUntilSolved(): Unit = {
      minMarkZeros()
      while (!solved) {
        relaxByMinUnmarked()
        resetMarks()
        minMarkZeros()
      }
    }

    private[onmf] class SolutionChooser() {
      val workerZeroJobsUnassigned: DenseVector[mutable.BitSet] = newWorkerZeroJobs()
      val jobZeroWorkersUnassigned: DenseVector[mutable.BitSet] = newJobZeroWorkers()
      val workerComparator: IntComparator = new SolutionChooser.AxisComparator(workerZeroJobsUnassigned)
      val jobComparator: IntComparator = new SolutionChooser.AxisComparator(jobZeroWorkersUnassigned)
      val workerIndices: Array[Int] = (0 until numWorkers).toArray
      val jobIndices: Array[Int] = (0 until numJobs).toArray
      val workerMarked: mutable.BitSet = mutable.BitSet((0 until numWorkers): _*) &~ workerUnmarked
      val jobMarked: mutable.BitSet = mutable.BitSet((0 until numJobs): _*) &~ jobUnmarked
      val workerAssignment: mutable.Map[Int, Int] =
        new mutable.HashMap[Int, Int](numWorkers, mutable.HashMap.defaultLoadFactor)
      val workerQ: IntHeapIndirectPriorityQueue =
        if (workerMarked.nonEmpty)
          new IntHeapIndirectPriorityQueue(workerIndices, workerMarked.toArray, numWorkers, workerComparator)
        else
          new IntHeapIndirectPriorityQueue(workerIndices, numWorkers, workerComparator)
      val jobQ: IntHeapIndirectPriorityQueue =
        if (jobMarked.nonEmpty)
          new IntHeapIndirectPriorityQueue(jobIndices, jobMarked.toArray, numJobs, jobComparator)
        else
          new IntHeapIndirectPriorityQueue(jobIndices, numJobs, jobComparator)

      def workerQFirstPriority: Int = {
        if (workerQ.isEmpty)
          Integer.MAX_VALUE
        else
          workerZeroJobsUnassigned(workerQ.first()).size
      }

      def jobQFirstPriority: Int = {
        if (jobQ.isEmpty)
          Integer.MAX_VALUE
        else
          jobZeroWorkersUnassigned(jobQ.first()).size
      }

      /**
       * Read off the minimum-cost solution (or one of such solutions).
       * Proceed one by one assigning marked workers and jobs to the first available single-marked 0 in their row or
       * column.
       * Priority order:
       * - For each marked worker or job, keep track of how many single-marked 0s in that row or column.
       * - Exclude any 0s double-marked by a crossing mark on any other column or row.
       * - Prioritize by count of single-marked 0s remaining for any worker or job, lowest count first.
       * - When workers and jobs tie for lowest count, choose a worker with that lowest count.
       * - As you assign a worker or job, remove it from any remaining sets of single-marked 0s for other workers or jobs.
       * - To remove the assigned worker or job, use the previously stored location of those zeros; no need to search
       * entire rows or columns.
       *
       * We use 2 indexed priority queues to store and prioritize the marked workers and jobs:
       * In particular, it.unimi.dsi.fastutil.ints.IntHeapIndirectPriorityQueue
       * "Indirect" priority queues give us operations beyond a normal heap:
       * - An index to access each item on the queue directly (not only the first item)
       * - IPQ.changed(index) method, to make the queue re-level any changed item in the heap.
       *
       * For their comparator, the priority queues use a worker's or job's count of remaining single-marked 0s.
       * As we complete each assignment, we remove the assigned worker or job from any others' bit-sets of single-marked
       * 0s, and call the change() method for its queue.
       */
      private[onmf] def chooseSolution(): Seq[Int] = {
        while (!workerQ.isEmpty || !jobQ.isEmpty) {
          if (workerQFirstPriority <= jobQFirstPriority) {
            val thisWorker = workerQ.dequeue()
            val jobAvailable = workerZeroJobsUnassigned(thisWorker) &~ jobMarked
            val job = jobAvailable.firstKey
            workerAssignment.put(thisWorker, job)
            workerMarked -= thisWorker
            val otherWorkersWithJob = workerMarked & jobZeroWorkersUnassigned(job)
            for (otherWorker <- otherWorkersWithJob) {
              workerZeroJobsUnassigned(otherWorker) -= job
              workerQ.changed(otherWorker)
            }
          } else {
            val thisJob = jobQ.dequeue()
            val workerAvailable = jobZeroWorkersUnassigned(thisJob) &~ workerMarked
            val worker = workerAvailable.firstKey
            workerAssignment.put(worker, thisJob)
            jobMarked -= thisJob
            val otherJobsWithWorker = jobMarked & workerZeroJobsUnassigned(worker)
            for (otherJob <- otherJobsWithWorker) {
              jobZeroWorkersUnassigned(otherJob) -= worker
              jobQ.changed(otherJob)
            }
          }
        }
        (0 until numWorkers).map(worker => workerAssignment(worker))
      }
    }

    private[onmf] object SolutionChooser {
      def apply(): SolutionChooser = {
        if (zeroMarking.solved) {
          new SolutionChooser()
        } else {
          zeroMarking.relaxUntilSolved()
          new SolutionChooser()
        }
      }

      private[onmf] class AxisComparator(axisZerosUnassigned: DenseVector[mutable.BitSet]) extends IntComparator {
        @Override
        def compare(a: Int, b: Int): Int = {
          Integer.compare(axisZerosUnassigned(a).size, axisZerosUnassigned(b).size)
        }
      }

    }
    def solutionChooser(): SolutionChooser = SolutionChooser()

  }
  private[onmf] object ZeroMarking {

  }

  val zeroMarking: ZeroMarking = ZeroMarking()

  def solve(): Unit = zeroMarking.relaxUntilSolved()

  override def toString: String = zeroMarking.toString
}

  object HungarianT {
    def apply(cost: DenseMatrix[Double]): HungarianT = {
      new HungarianT(cost)
    }

    /**
     *
     */
    def main(args: Array[String]): Unit = {
//      val m4by4 = DenseMatrix((1 to 16).map(_.toDouble).toArray).reshape(4, 4).t
//      val m4by4sqr = m4by4 *:* m4by4
//      var mExpand = m4by4sqr.copy
//      mExpand = DenseMatrix.horzcat(mExpand(0 until mExpand.rows, 0 to 0), mExpand)
//      mExpand = DenseMatrix.vertcat(mExpand((0 to 0), (0 until mExpand.cols)), mExpand)
//      //      mExpand = mExpand(Seq(0) ++ (0 until mExpand.rows), (0 until mExpand.cols)).toDenseMatrix
//      //      mExpand = mExpand((0 until mExpand.rows), Seq(0) ++ (0 until mExpand.cols)).toDenseMatrix
//      val h = HungarianT(cost = mExpand)
//      println(h.state.toString)
//      h.solve()
//      val A = h.assignment()
//      A.assign()
//      println(A)
//      val mtri0 = DenseMatrix((0.0, 0.0, 1.0, 1.0, 1.0),
//        (0.0, 0.0, 0.0, 1.0, 1.0),
//        (1.0, 0.0, 0.0, 0.0, 1.0),
//        (1.0, 1.0, 0.0, 0.0, 0.0),
//        (1.0, 1.0, 1.0, 0.0, 0.0))
//      val h2 = HungarianT(cost = mtri0)
//      println(h2.state.toString)
//      h2.solve()
//      val a2 = h2.assignment()
//      a2.assign()
      val mtri0tall = DenseMatrix((0.0, 0.0, 1.0, 1.0, 1.0),
        (0.0, 0.0, 0.0, 1.0, 1.0),
        (1.0, 0.0, 0.0, 0.0, 1.0),
        (1.0, 1.0, 0.0, 0.0, 0.0),
        (1.0, 1.0, 1.0, 0.0, 0.0),
        (1.0, 1.0, 1.0, 0.0, 0.0),
      )
      val h3 = HungarianT(cost = mtri0tall)
      println(h3.zeroMarking.toString)
      h3.zeroMarking.relaxUntilSolved()

    }
  }