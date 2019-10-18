package com.deanelzinga.kuhnmunkres {

  import breeze.linalg._
  import com.deanelzinga.kuhnmunkres.Hungarian._
  import it.unimi.dsi.fastutil.ints.{IntComparator, IntComparators, IntHeapIndirectPriorityQueue, IntHeapSemiIndirectPriorityQueue}
  import org.roaringbitmap.RoaringBitmap
  import scala.collection.{immutable, mutable}
  import java.util.Comparator

  import scala.collection.mutable.HashMap

  /**
   * Kuhn-Munkres typically stated with workers in rows and jobs in columns. Consider
   * transposing this for Breeze, since columns are dominant. Requires reorienting
   * workBoard and the order, (::, *) vs (::, *), of indexing matrix operations.
   * *
   * Copy cost matrix to mutable, markable cost slate:
   *- Transpose cost matrix as needed: shorter axis indexes worker; longer axis, jobs.
   *- Prepare sets of unmarked workers and jobs, so we can mark columns or rows containing 0.
   *- We subtract each worker's minimum job-cost from all the job costs for that worker.
   *- We subtract each job's minimum worker cost from all the worker costs for that job.
   * STEP 1. For each worker, SUBTRACT the min job cost from all the job costs for that worker.
   * For each worker, subtract their min job cost from all their job costs:
   * STEP 2. For each job, SUBTRACT the min worker cost from all the worker costs for that job.
   * (Note, for jobs with a worker cost set to 0 in Step 1, subtracting 0 obviously has no effect.)
   * For each job, subtract its min worker costs from all their worker costs:
   * *
   * STEP 3. Starting from an unmarked, transformed cost matrix...
   *- Mark all the 0s of the transformed cost matrix with a minimum number of marks on
   * workers or jobs (rows or columns).
   * Find and make next mark, if any:
   * At each step, we use a greedy algorithm:
   *- Find the max count of unmarked, 0-cost jobs for any unmarked worker.
   *- Find the maximum number of unmarked, 0-cost workers for any unmarked job.
   *- If there are no more unmarked, 0-cost jobs, go to STEP 4.
   *- If either is greater, choose it. If they are equal, prefer the worker.
   *- Mark that worker or job (remove it from unmarked workers or jobs).
   *- Remove that worker or job from all jobs' or workers' sets of 0-cost workers or jobs.
   * *
   * STEP 4. TEST FOR OPTIMALITY: (i) If the minimum of 0-covering marks is N,
   * the number of workers, then an optimal assignment of 0s is possible and "we are
   * finished"-- Go to Step 6.
   * *
   * STEP 5. Find the minimum cost entry still unmarked.
   *- Subtract this entry from each unmarked worker.
   *- Add this entry to each MARKED job. Return to STEP 3.
   * *
   * STEP 6: Read off an available assignment from the covered 0s (not completely trivial).
   */
  class Hungarian(cost: DenseMatrix[Double]) {

    // Originally the algorithm was stated with rows fewer than columns. This matters for some of the order
    // of operations during the algorithm. For Breeze, which is column oriented, it probably makes more sense
    // for columns to be longer than rows. In anticipation of this change, the shorter axis is called workers
    // and the longer axis is called jobs, to separate the algorithmic decisionmaking from the implementation
    // details.
    protected var costT: DenseMatrix[Double] = {
      if (cost.rows <= cost.cols)
        cost.toDenseMatrix
      else
        cost.t.toDenseMatrix
    }


    protected var costX: DenseMatrix[Double] = costT.copy
    {
      val jobMinWorkerCost = min(costX(*, ::))
      costX(::, *) :-= jobMinWorkerCost
    }

    {
      val workerMinJobCost = min(costX(::, *)).t
      costX(*, ::) :-= workerMinJobCost
    }


    def numWorkers: Int = costX.rows
    def numJobs: Int = costX.cols
    private def newWorkerUnmarked() = mutable.BitSet(0 until numWorkers: _*)
    private def newJobUnmarked() = mutable.BitSet(0 until numJobs: _*)

    private def newWorkerZeroJobs(): DenseVector[mutable.BitSet] = costX(*, ::).
      map(c => where(c :== 0.0)).
      map(mutable.BitSet(_: _*))

    private def newJobZeroWorkers(): DenseVector[mutable.BitSet] = costX(::, *).
      map(c => where(c :== 0.0)).t.
      map(mutable.BitSet(_: _*))

    case class Mark(worker: Boolean, index: Int)

    /**
     *
     * @param costX transformned cost matrix, with workers the shortest axis
     * @param workerUnmarked
     * @param jobUnmarked
     * @param workerZeroJobs
     * @param jobZeroWorkers
     */
    case class State(costX: DenseMatrix[Double] = costX,
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
      def reset(): Unit = {
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
        val jobIndices = (0 until numJobs) map (job => f"$job%4d") mkString("    " + " job", "", "\n")
        val jobZeroSizes = (0 until numJobs) map (
          job => {
            val numZeroWorkers = jobZeroWorkers(job).size
            if (jobUnmarked(job)) "   |" else f"$numZeroWorkers%4d"
          }) mkString(" wrk" + " mrk", "", "\n")
        val lines = (0 until numWorkers) map {
          worker => {
            val numZeroJobs = workerZeroJobs(worker).size
            Seq(f"$worker%4d",
              if (workerUnmarked(worker)) "  --" else f"$numZeroJobs%4d",
              costX(worker, ::).t.toScalaVector().map(c => if (c == 0.0) "   ." else f"$c%4.0f").mkString
            ).mkString
          }
        }
        Seq(jobIndices, jobZeroSizes, lines.mkString("", "\n", "\n")).mkString
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

      def markAllZeros(): Unit = {
        var markMaybe = nextMarkOption
        while (markMaybe.nonEmpty) {
          state.markZeroIfAny(markMaybe)
          markMaybe = state.nextMarkOption
        }
      }

      /** number of workers > Number of marks */
      def solved: Boolean = {
        numWorkers == (numWorkers - workerUnmarked.size) + (numJobs - jobUnmarked.size)
      }

      def reduceByMinUnmarked(): Unit = {
        if (workerUnmarked.nonEmpty && jobUnmarked.nonEmpty) {
          val minUnmarked = min(costX(workerUnmarked.toSeq, jobUnmarked.toSeq))

          // Subtract min unmarked cost from every unmarked worker;
          costX(workerUnmarked.toSeq, ::) :-= minUnmarked

          // Add min unmarked cost to every MARKED column:
          costX(::, (0 until numJobs).filterNot(jobUnmarked(_))) :+= minUnmarked
        }
        else {}
      }

      // FIXME: Refactor Hungarian.solve() to be Hungarian.state.solve()
      def solve(): Unit = {
        while(!solved) {
          reduceByMinUnmarked()
          reset()
          markAllZeros()
        }
      }
    }

    object State {

    }

    def solve(): Unit = {
      while (!state.solved) {
        state.reduceByMinUnmarked()
        state.reset() // Reset all marks.
        state.markAllZeros()
      }
    }

    var state: State = State()
    state.markAllZeros()

    private[kuhnmunkres] class Assignment(state: State) {

      /**
       *
       *
       * @param axisZerosUnassigned
       */
      class AxisComparator(axisZerosUnassigned: DenseVector[mutable.BitSet]) extends IntComparator {
        @Override
        def compare(a: Int, b: Int): Int = {
          Integer.compare(axisZerosUnassigned(a).size, axisZerosUnassigned(b).size)
        }
      }

      // protected var numUnexpectedStates: Long = 0L
      val workerZeroJobsUnassigned: DenseVector[mutable.BitSet] = newWorkerZeroJobs()
      val jobZeroWorkersUnassigned: DenseVector[mutable.BitSet] = newJobZeroWorkers()
      val workerComparator: IntComparator = new AxisComparator(workerZeroJobsUnassigned)
      val jobComparator: IntComparator = new AxisComparator(jobZeroWorkersUnassigned)
      val workerIndices: Array[Int] = (0 until numWorkers).toArray
      val jobIndices: Array[Int] = (0 until numJobs).toArray
      val workerMarked: mutable.BitSet = mutable.BitSet((0 until numWorkers): _*) &~ state.workerUnmarked
      val jobMarked: mutable.BitSet = mutable.BitSet((0 until numJobs): _*) &~ state.jobUnmarked
      // val workerAssigned: mutable.BitSet = new mutable.BitSet(initSize = numWorkers)  // Consider removing
      // val jobAssigned: mutable.BitSet = new mutable.BitSet(initSize = numJobs) // Consider removing
      val workerAssignment: mutable.Map[Int, Int] =
        new mutable.HashMap[Int, Int](numWorkers, mutable.HashMap.defaultLoadFactor)
      val jobAssignment: mutable.Map[Int, Int] =
        new mutable.HashMap[Int, Int](numWorkers, mutable.HashMap.defaultLoadFactor) // short axis length
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
      def anyQNonEmpty: Boolean = {
        !workerQ.isEmpty || !jobQ.isEmpty
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
      private[kuhnmunkres] def assign(): Unit = {
        while (anyQNonEmpty) {
          if (workerQFirstPriority <= jobQFirstPriority) {
            val thisWorker = workerQ.dequeue()
            // workerAssigned += thisWorker
            val jobAvailable = workerZeroJobsUnassigned(thisWorker) &~ jobMarked
            val job = jobAvailable.firstKey
            // jobAssigned += job
            workerAssignment.put(thisWorker, job)
            jobAssignment.put(job, thisWorker)
            workerMarked -= thisWorker
            val otherWorkersWithJob = workerMarked & jobZeroWorkersUnassigned(job)
            for (otherWorker <- otherWorkersWithJob) {
              workerZeroJobsUnassigned(otherWorker) -= job
              workerQ.changed(otherWorker)
            }
          } else {
            val thisJob = jobQ.dequeue()
            // jobAssigned += thisJob
            val workerAvailable = jobZeroWorkersUnassigned(thisJob) &~ workerMarked
            val worker = workerAvailable.firstKey
            // workerAssigned += worker
            jobAssignment.put(thisJob, worker)
            workerAssignment.put(worker, thisJob)
            jobMarked -= thisJob
            val otherJobsWithWorker = jobMarked & workerZeroJobsUnassigned(worker)
            for (otherJob <- otherJobsWithWorker) {
              jobZeroWorkersUnassigned(otherJob) -= worker
              jobQ.changed(otherJob)
            }
          }
        }
      }
    }

    private[kuhnmunkres] object Assignment {
      def apply(): Assignment = {
        if (state.solved) {
          new Assignment(state)
        } else {
          state.solve()
          new Assignment(state)
        }
      }
    }

    def assignment(): Assignment = Assignment()
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
//      val h = Hungarian(cost = mExpand)
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
//      val h2 = Hungarian(cost = mtri0)
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
      val h3 = Hungarian(cost = mtri0tall)
      println(h3.state.toString)
      h3.state.solve()
      val a3 = h3.assignment()
      a3.assign()

    }
  }
}