package com.deanelzinga.onmf {

  import breeze.linalg._
  import it.unimi.dsi.fastutil.ints.{IntComparator, IntHeapIndirectPriorityQueue}

  import scala.collection.mutable

  /**
   * Kuhn-Munkres typically stated with workers in rows and jobs in columns. Consider
   * transposing this for Breeze, since columns are dominant. Requires reorienting
   * workBoard and the order, (::, *) vs (::, *), of indexing matrix operations.
   * *
   * Copy cost matrix to mutable, markable cost matrix:
   *- Transpose cost matrix as needed: call shorter axis, worker; longer axis, jobs.
   *- Prepare sets of unmarked workers and jobs, so we can mark columns or rows containing 0.
   *- We subtract each worker's minimum job-cost from all the job costs for that worker.
   *- We subtract each job's minimum worker cost from all the worker costs for that job.
   * STEP 1. REDUCE EACH WORKER: For each worker, SUBTRACT the min job cost from all the job costs for that worker.
   * For each worker, subtract their min job cost from all their job costs:
   * STEP 2. REDUCE EACH JOB: For each job, SUBTRACT the min worker cost from all the worker costs for that job.
   * (Note, for jobs with a worker cost set to 0 in Step 1, subtracting 0 obviously has no effect.)
   * For each job, subtract its min worker costs from all their worker costs:
   * *
   * STEP 3. MARK ALL ZEROS: Starting from an unmarked, transformed cost matrix,
   * mark all the 0s of the transformed cost matrix with the minimum number of marks on
   * workers or jobs (rows or columns).
   *
   * While any 0 cells remain unmarked:
   * 3A: Find and mark the next mark:
   * - Over all unmarked workers, find the max count per worker of unmarked, 0-cost cells (jobs).
   * - Over all unmarked jobs, find the max count per job of unmarked, 0-cost cells (workers).
   * If either of the max counts is 0, then they both are. We are done marking. Go to step 4.
   * Otherwise both max counts are greater than 0:
   * If either of the max counts (per worker or per job) is strictly greater than the other, choose that axis.
   * Mark the first worker or job with that greater max count. Remove the worker or job from unmarked jobs' or workers'
   * list of 0 cells. Go on to the next mark (3A).
   * Otherwise the max counts are equal. For workers and jobs, count the number of each having ANY
   * unmarked 0 cells.
   * - If either workers or jobs has a strictly smaller count, then that axis is the shorter route to completion.
   * Choose the axis with the smaller count. Mark the first line of cells on that axis having the max count of
   * unmarked 0s. Remove the marked worker or job from unmarked jobs' or workers' list of 0 cells. Go on to the next mark (3A).
   * - Otherwise the counts are equal. In this case, choose the shorter axis (workers). Mark the first worker having
   * the maximum available unmarked 0 cells. Remove the marked worker or job from unmarked jobs' or workers' list of 0 cells.
   * Go on to the next mark (3A).
   *
   * STEP 4. TEST FOR OPTIMALITY: (i) If the minimum of 0-covering marks is N,
   * the number of workers, then an optimal assignment of 0s is possible and "we are
   * finished"-- Go to Step 6.
   * Otherwise the minimum of 0-covering marks is less than N. Go to STEP 5.
   * *
   * STEP 5. RELAXATION: Find the minimum cost entry still unmarked.
   *- Subtract this entry from each unmarked worker.
   *- Add this entry to each MARKED job. Return to STEP 3.
   * *
   * STEP 6: Read off an available assignment from the covered 0s (not completely trivial).
   */
  class KuhnMunkres(cost: DenseMatrix[Double]) {
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
    def reduceRows(): Unit =
    {  // Reduce rows
      val jobMinWorkerCost = min(costX(*, ::))
      costX(::, *) :-= jobMinWorkerCost
    }
    reduceRows()
    def reduceCols(): Unit = {  // Reduce columns
      val workerMinJobCost = min(costX(::, *)).t
      costX(*, ::) :-= workerMinJobCost
    }
    reduceCols()
    def numWorkers: Int = costX.rows
    def numJobs: Int = costX.cols
    private def newWorkerUnmarked() = mutable.BitSet(0 until numWorkers: _*)
    private def newJobUnmarked() = mutable.BitSet(0 until numJobs: _*)

    private def newWorkerZeroJobs(): DenseVector[mutable.BitSet] =  { costX(*, ::).
      map(c => where(c :== 0.0)).
      map(mutable.BitSet(_: _*))
    }
    private def newJobZeroWorkers(): DenseVector[mutable.BitSet] = { costX(::, *).
      map(c => where(c :== 0.0)).t.
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
          val minUnmarked = min(costX(workerUnmarked.toSeq, jobUnmarked.toSeq))

          // Subtract min unmarked cost from every unmarked worker;
          costX(workerUnmarked.toSeq, ::) :-= minUnmarked

          // Add min unmarked cost to every MARKED column:
          costX(::, (0 until numJobs).filterNot(jobUnmarked(_))) :+= minUnmarked
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

  object KuhnMunkres {
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
}