package com.deanelzinga.onmf {

  import breeze.linalg._
  import breeze.numerics._
  import it.unimi.dsi.util._

  import scala.util.Random
  import Gopa._

  /**
   * Greedy Orthogonal Pivoting Algorithm for Non-Negative Matrix Factorization (GOPA-NMF).
   * http://proceedings.mlr.press/v97/zhang19r.html
   * - Currently implements the sequential algorithm (Algorithm 1)
   * @param x (Matlab `X`) main observation matrix to factor into nonnegative W and H.
   * @param gnd (Matlab `gnd`) I think for "ground truth", since it's used to calculate Ac (accuracy?)
   * @param k (Matlab `k`) dimensions to reduce to.
   * @param numT (Matlab `T`) number of full iterations.
   * @param ratio (Matlab `ratio`)
   */
  class Gopa(x: DenseMatrix[Double], gnd: DenseVector[Int], k: Int, numT: Int, ratio: Double) {
    val n: Int = x.rows
    var w0: DenseMatrix[Double] = randOrtho(x.rows, k)
    // var ac: DenseVector[Double] = DenseVector.zeros[Double](numT) // FIXME
    var count = 1
    for (t <- 0 until numT) {
      println("t: " + t)

      // Cluster label: column index of the nonzero element (the max, since the rest are 0): (argmax(w0, Axis._1))
      val res: DenseVector[Int] = argmax(w0, Axis._1)  // Cluster assignments [~,res] = max(W0');
      // ac(t) = Gopa.accuracy(gnd, res)  // FIXME
      val h0: DenseMatrix[Double] = w0.t * x       // Update matrix H by H=W′X;
      count = count + 1
      val r: DenseMatrix[Double] = x * h0.t        // R = X*H0';
      normalizeCol(w0)  // Rescale W columns p and q to norm 1; W0 = W0*diag(1./sqrt(1e-10+sum(W0.*W0)));
      val u: Seq[Int] = randPerm(n)
      for (j <- 0 until n) {
        println("j: " + j)
        if (Gopa.rand() <= ratio) {
          val v = gop1(w0, r, u(j))
          w0 = v.w0
        }
      }
    }

    case class Gop1(w0: DenseMatrix[Double] = w0,
                    e: Int = 0,        // e = q;..e = find(T == mini);..e = e(1);..if(e~=q);..
                    status: Int = -1,  // %status  % -1: error; % 0: no potential; % 1: can reduce
                    o: Double = 0.0,   // [o,Wpn,Wqn,x] = sml(W0,R,q,e,l0);..[o,Wqn,x] = sml_0(W0,R,q,l0);.. o = W0(l,q);
                    t: DenseVector[Double] = DenseVector.zeros[Double](k),   // T = zeros(1,k);...T(i) = obj;
                    p: DenseVector[Double] = DenseVector.zeros[Double](k))   // P = zeros(1,k);...P(i) = x;

    case class Sml(dif: Double = 0.0,
                   wpn: DenseVector[Double] = DenseVector.zeros[Double](n),
                   wqn: DenseVector[Double] = DenseVector.zeros[Double](n),
                   x: Double = .5)

    case class Sml_0(dif: Double = 0.0,
                     wqn: DenseVector[Double] = DenseVector.zeros[Double](n),
                     x: Double = .5)

    def gop1(w0: DenseMatrix[Double], r: DenseMatrix[Double], l0: Int): Gop1 = {
      val k: Int = w0.cols
      val w: DenseVector[Double] = w0(l0, ::).t
      val qs: BitVector = (w >:> 0.0)
      if (qs.activeSize == 0) {
        Gop1(w0, e = 0, status = -1, o = 0.0, t = randVector(k), p = DenseVector.zeros[Double](k))
      } else {
        val q = qs.activeKeysIterator.next()
        val t = DenseVector.zeros[Double](k)
        val p = DenseVector.zeros[Double](k)

        // Swap the non-zero entry of W[i,:] from the original location (qth column) to each of the {1,2,...,k} locations;
        for (i <- 0 until k) {  // For
          //
          val v =
            if (i != q)
              sml(w0, r, q, i, l0)
            else
              sml_0(w0, r, q, l0)
          t(i) = v.dif  // Store the objective change for each candidate column.
          p(i) = v.x    // Store the optimal x for each candidate column.
        }
        val dex = (t <:< 0.0) &:& (p >:> 0.0) // Indices where objective improves and optimal switch is positive
        if (dex.activeSize == 0) {  // %print('total zero row, warning');
          Gop1(w0 = w0, e = q, status = 0, o = 0.0, t, p)
        } else {
          val mini = min(t(dex))              // Best objective improvement
          val e = (t :== mini).activeKeysIterator.next()  // Index of best objective improvement
          if (e != q ) {  // 2-column case: swap q to different column, e
            val v = sml(w0, r, q, e, l0)
            w0(::, e) := v.wpn
            w0(::, q) := v.wqn
            Gop1(w0, e , status = 1, o = v.dif, t, p)
          } else {
            val v: Sml = sml_0(w0, r, q, l0)
            w0(::, q) := v.wqn
            Gop1(w0, e, status = 1, o = v.dif, t, p)
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
        if (index < cols) // Ensure every column has at least 1 non-zero value (1.0)
          m(randRow(index), index) = 1.0 // Sequential instead of random for the first few rows. We will rearrange later.
        else
          m(randRow(index), random.nextInt(cols)) = 1.0
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