package com.deanelzinga.onmf {

  import breeze.linalg._
  import breeze.numerics._
  import com.deanelzinga.onmf.Gopa._
  import it.unimi.dsi.util._

  import scala.util.Random

  /**
   * Greedy Orthogonal Pivoting Algorithm for Non-Negative Matrix Factorization (GOPA-NMF).
   * http://proceedings.mlr.press/v97/zhang19r.html
   * - Currently implements the sequential algorithm (Algorithm 1)
   * @param x (Matlab `X`) main observation matrix to factor into nonnegative W and H.
   * @param gnd (Matlab `gnd`) I think for "ground truth", since it's used to calculate Ac (accuracy?)
   * @param k (Matlab `k`) dimensions to reduce to.
   * @param numT (Matlab `T`) number of full iterations.
   * @param ratio (Matlab `ratio`) ratio of whether to update or leave unchanged each row of each iteration.
   */
  class Gopa(x: DenseMatrix[Double], gnd: DenseVector[Int], k: Int, numT: Int, ratio: Double) {
    // Initialize: TODO: Consider moving to method init()
    val n: Int = x.rows
    var w0: DenseMatrix[Double] = randOrtho(x.rows, k)

    // Optimization, main loop. TODO: Consider moving to method fit() or optimize()
    // var ac: DenseVector[Double] = DenseVector.zeros[Double](numT) // FIXME: Reverse engineer ac, res use case.
    var count = 1
    for (t <- 0 until numT) {
      println("t: " + t)

      // Cluster label: column index of the nonzero element (the max, since the rest are 0): (argmax(w0, Axis._1))
      // val res: DenseVector[Int] = argmax(w0, Axis._1)  // Cluster assignments [~,res] = max(W0'); // FIXME: res, ac
      // ac(t) = Gopa.accuracy(gnd, res)  // FIXME: figure out what ORNF
      val h0: DenseMatrix[Double] = w0.t * x       // Update matrix H by H=W′X;
      count = count + 1
      val r: DenseMatrix[Double] = x * h0.t        // R = X*H0';
      val u: Seq[Int] = randPerm(n)
      for (j <- 0 until n) {
        println("j: " + j)
        if (Gopa.rand() <= ratio) {
          gop1(w0, r, u(j))  // FIXME: Change to mutate operation and Unit return value.
        }
      }
    }

    /**
     *
     * @param w0
     * @param p
     * @param status
     * @param obj objective change of swapping from column p to q (or leaving at p), rescaling to optimal value.
     * @param objs
     * @param x1s
     */
    case class Gop1(w0: DenseMatrix[Double] = w0,
                    p: Int = 0, // e = q;..e = find(T == mini);..e = e(1);..if(e~=q);..
                    status: Int = -1, // %status  % -1: error; % 0: no potential; % 1: can reduce
                    obj: Double = 0.0, // [o,Wpn,Wqn,x] = sml(W0,R,q,e,l0);..[o,Wqn,x] = sml_0(W0,R,q,l0);.. o = W0(l,q);
                    objs: DenseVector[Double] = DenseVector.zeros[Double](k), // T = zeros(1,k);...T(i) = obj;
                    x1s: DenseVector[Double] = DenseVector.zeros[Double](k))   // P = zeros(1,k);...P(i) = x;

    case class Sml(dif: Double = 0.0,
                   wpn: DenseVector[Double] = DenseVector.zeros[Double](n),
                   wqn: DenseVector[Double] = DenseVector.zeros[Double](n),
                   x1: Double = .5)

    case class Sml_0(dif: Double = 0.0,
                     wqn: DenseVector[Double] = DenseVector.zeros[Double](n),
                     x1: Double = .5)

    /**
     *
     * @param w0
     * @param r
     * @param l0
     * @return
     */
      // FIXME: Change to Unit type and change w0 in place.
      // FIXME: Handle renormalize here only for column p; and other column q, if any.
    def gop1(w0: DenseMatrix[Double], r: DenseMatrix[Double], l0: Int): Gop1 = {
      val k: Int = w0.cols
      val w: DenseVector[Double] = w0(l0, ::).t
      val qs: BitVector = (w >:> 0.0)
      if (qs.activeSize == 0) {
        () // Gop1(w0, p = 0, status = -1, obj = 0.0, objs = randVector(k), x1s = DenseVector.zeros[Double](k))
      } else {
        val q = qs.activeKeysIterator.next()  // Get first integer in qs (indices of 1s), considered as a set
        val objs = DenseVector.zeros[Double](k)  // Optimal objective improvement for each column i of [0, k-1]
        val x1s = DenseVector.zeros[Double](k)  // Optimal value to change w0(l0,x1s) to when swapping to optimal location.

        // Swap the non-zero entry of W[i,:] from the original location (qth column) to each of the {1,2,...,k} locations;
        // TODO: Save the optimal i and vectors wpn and wqn
        for (p <- 0 until k) {  // For
          // v is struct with full matrix, optimal column, status,
          val v =
            if (p != q)
              sml(w0, r, q, p, l0)
            else
              sml_0(w0, r, q, l0)
          objs(p) = v.dif  // Store the objective change for each candidate column.
          x1s(p) = v.x1    // Store the optimal x for each candidate column.
        }
        val pPossible = (objs <:< 0.0) &:& (x1s >:> 0.0) // Indices where objective improves and optimal switch is positive
        if (pPossible.activeSize == 0) {  // %print('total zero row, warning');
          ()  // Gop1(w0 = w0, p = q, status = 0, obj = 0.0, objs, x1s)
        } else {
          // TODO: Here we revisit the optimal index and recompute the vectors, when we could have saved them.
          val minObj = min(objs(pPossible))              // Best objective improvement
          val p = (objs :== minObj).activeKeysIterator.next()  // First index equalling best objective
          if (p != q ) {  // 2-column case: swap q to different column, p
            val v = sml(w0, r, q, p, l0)
            w0(::, p) := v.wpn
            w0(::, q) := v.wqn
            normalizeCol(w0, p)
            normalizeCol(w0, q)
            // Gop1(w0, p , status = 1, obj = v.dif, objs, x1s)
          } else {
            val v: Sml = sml_0(w0, r, q, l0)
            w0(::, q) := v.wqn
            normalizeCol(w0, p)
            // Gop1(w0, p, status = 1, obj = v.dif, objs, x1s)
          }
        }
      }
    }

    /**
     * From function [dif,Wpn,Wqn,x]= sml(W0,R,q,p,l);
     * @param w0
     * @param r
     * @param q
     * @param p
     * @param l
     * @return
     */
    def sml(w0: DenseMatrix[Double], // w0: left factor so far, cluster assignments
            r: DenseMatrix[Double],  // R = X*H0'; precomputed dot products of data (X) and cluster-centers (H0)
            q: Int,  // Column q: unique non-zero value in row l.
            p: Int,  // Column p: candidate destination for swap, different from q
            l: Int  // Row l (elsewhere l0):
           ): Sml = {
      val o: Double = w0(l, q)                 // In row l, o is that non-zero value of w0 at (l, q)
      val wp: DenseVector[Double] = w0(::, p)  // wp: p candidate column of w0
      val rp: DenseVector[Double] = r(::, p)   // rp: p candidate column of precomputed dot products of X with cluster centers
      val ip: BitVector = wp >:> 0.0             // ip: indices where wp > 0
      val wcp: DenseVector[Double] = wp(ip).copy  // wcp: nonzero slice of candidate destination column

      val wq: DenseVector[Double] = w0(::, q)  // wq: q current non-zero column of w0
      val rq: DenseVector[Double] = r(::, q)   // rq: q current non-zero column of r

      val a: Double = wcp dot rp(ip) // = wcp.t * rp(ip); rp(ip) is rp restricted to indices where w0 > 0.
      val b: Double = rp(l)          // Row l of the column p of r
      val x: Double = b / sqrt(a*a + b*b)

      val y: Double = sqrt(1.0 - x*x)  // Amount to mul
      val wpn: DenseVector[Double] = wp.copy
      val wqn: DenseVector[Double] = wq.copy
      wpn(l) = x
      wpn(ip) := wcp *:* y
      wqn(l) = 0.0
      wqn :/= sqrt(1.0 - o*o)  // Divide in place;..Wqn = Wqn/sqrt(1-o^2)
      val dif = -(a*(y - 1.0) + x*rp(l) + ((wqn - wq) dot rq))  // dif = -(a*(y-1) + x*Rp(l) + (Wqn-Wq)'*Rq);
      Sml(dif, wpn, wqn, x)
    }

    /**
     * Same-column case, column q. Calculates the change in objective (potential improvement) from a same-column rescale
     * of q(l) and renormalization of column q.
     *
     * From sml_0.m: function [dif,Wqn,x]= sml_0(W0,R,q,l);
     * - Altered to align the signature of the return values (4-tuple) with the 2-column case, sml(), to align the
     * signaturassignments of the function return values.
     * @param w0
     * @param r
     * @param q
     * @param l
     * @return
     */
    def sml_0(w0: DenseMatrix[Double],
              r: DenseMatrix[Double],
              q: Int,
              l: Int
             ): Sml = {
      val wq: DenseVector[Double] = w0(::, q)     // Wq = W0(:,q);
      val rq: DenseVector[Double] = r(::, q)      // Rq = R(:,q);
      val o: Double = wq(l)                       // o = Wq(l);
      val wqn: DenseVector[Double] = wq.copy      // Wqn = Wq;

      wq(l) = 0.0                                 // Wq(l) = 0;
      val iq: BitVector = (wq >:> 0.0)            // I = find(Wq>0);
      val wcq: SliceVector[Int, Double] = wq(iq)  // Wcq = Wq(I);
      val rcq: SliceVector[Int, Double] = rq(iq)  // Rcq = Rq(I);
      val b: Double = rq(l)                       // b = Rq(l);

      val j0: Double = -(wcq dot rcq + o*b)       // J0 = -(Wcq'*Rcq + o*b);

      val wcqp: DenseVector[Double] = wcq / sqrt(1.0 - o*o)
      val a: Double = wcqp dot rcq
      val x: Double = b / sqrt(a*a + b*b)
      val y: Double = sqrt(1.0 - x*x)
      val j1: Double = -(a*y + x*b)
      val dif: Double = j1 - j0

      wqn(l) = x
      wqn(iq) := wcqp * y
      Sml(dif, wqn, wqn, x)  // Same return types and arity as sml()
    }

  }

  object Gopa {
    var random: Random = new XoRoShiRo128PlusPlusRandom()
    def rand(): Double = random.nextDouble()

    /**
     * Normalizes in place all columns of a matrix.
     * @param m
     */
    def normalizeCol(m: DenseMatrix[Double]): Unit = {
      m := m * diag(1.0 /:/ sqrt(sum(m *:* m, Axis._0).t))
    }

    /**
     * Normalizes in place only column q of a matrix.
     * @param m
     * @param q
     */
    def normalizeCol(m: DenseMatrix[Double], q: Int): Unit = {
      val mq: DenseVector[Double] = m(::, q)
      m(::, q) :/= sqrt(mq dot mq)
    }
    def randPerm(n: Int): Seq[Int] = {
      random.shuffle((0 until n).toVector)
    }

    def randOrtho(rows: Int, cols: Int): DenseMatrix[Double] = {
      val m = DenseMatrix.zeros[Double](rows, cols)
      val randRow = random.shuffle((0 until rows).toVector)

      0 until rows foreach { index =>
        val x = 1.0
        // Gaussian(0.0, sqrt(cols/rows))  // TODO: Consider using this random value instead of 1 to make these
          // Typically norm of vector of N Gaussian(0,1) is ~sqrt(N). With sigma=sqrt(cols/rows), this makes each
          // column approx norm 1.0 even before normalization. Using Gaussian variates instead of 1.0 might give
          // these some of the representation power of Random Projections, a separate technique following from
          // Johnson-Lindenstrauss Lemma.
        if (index < cols) // Ensure every column has at least 1 non-zero value (1.0)
          m(randRow(index), index) = x
        else
          m(randRow(index), random.nextInt(cols)) = x
      }

      // Normalize columns:
      normalizeCol(m)
      m
    }

    def randVector(m: Int): DenseVector[Double] = {
      val v = DenseVector.zeros[Double](m)
      v(random.nextInt(m)) = random.nextDouble()
      v
    }

    def bestMap(gnd: DenseVector[Int], res: DenseVector[Int]): DenseVector[Int] = {
      res // FIXME: stub to get to compile, replace with KuhnMunkres of clustering assignments
    }

    def accuracy(gnd: DenseVector[Int], res: DenseVector[Int]): Double = {
      val bestRes = bestMap(gnd, res)
      (bestRes :== gnd).activeSize.toDouble / gnd.length
    }

    def main(args: Array[String]): Unit = {
      // val shuffled = random.shuffle((0 until 10).toVector)

      val randOrtho1: DenseMatrix[Double] = randOrtho(10, 4)
      val iWhere: Double = I(randOrtho1 == 0.0)
      println(randOrtho1)
      // shuffled.foreach( x => println(x) )
    }
  }

}