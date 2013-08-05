package com.nr.min;


import static com.nr.NRUtil.*;
import static java.lang.Math.*;
import com.nr.RealValueFun;

/**
 * Downhill Simplex Method in Multidimensions
 * 
 * downhill simplex minimization
 * Copyright (C) Numerical Recipes Software 1986-2007
 * Java translation Copyright (C) Huang Wen Hui 2012
 *
 * @author hwh
 *
 */
// Multidimensional minimization by the downhill simplex method of Nelder and Mead
public class Amoeba {
  final double ftol;
  int nfunc;   // The number of function evaluations.
  int mpts;  
  int ndim;
  double fmin;   // Function value at the minimum.
  double[] y;   // Function values at the vertices of the simplex.
  double[][] p;    // Current simplex
  
  // The constructor argument ftoll is the fractional convergence tolerance to be 
  // achieved in the function value
  public Amoeba(final double ftoll) {
    this.ftol = ftoll;
  }
  
  /*
   Multidimensional minimization of the function f(x), where x[0..ndim-1] 
   is a vector in ndim dimensions, by the downhill simplex method of Nedler and Mead.
   The initial simplex is specified as in equation (10.5.1) by a point[0..ndim-1] and a 
   constant displacement del along each coordinate direction. Returned is the location of the minimum.
   */
  public double[] minimize(final double[] point, final double del, final RealValueFun func) {
    double[] dels = buildVector(point.length,del);
    return minimize(point,dels,func);
  }
  
  // Alternative inteface that takes different displacements dels[0..ndim-1] in different directions
  // for the initial simplex.
  public double[] minimize(final double[] point, final double[] dels, final RealValueFun func) {
    int ndim=point.length;
    double[][] pp = new double[ndim+1][ndim];
    for (int i=0;i<ndim+1;i++) {
      for (int j=0;j<ndim;j++)
        pp[i][j]=point[j];
      if (i !=0 ) pp[i][i-1] += dels[i-1];
    }
    return minimize(pp,func);
  }

  // Most general interface: initial simplex specified by the matrix pp[0..ndim][0..ndim-1]
  // Its ndim+1 rows are ndim-dimensional vectors that are vertices of the starting simplex
  public double[] minimize(final double[][] pp, final RealValueFun func) {
    final int NMAX=5000;    // Maximum allowed number of function evaluations
    final double TINY=1.0e-10;
    int ihi,ilo,inhi;
    mpts=pp.length;
    ndim=pp[0].length;
    double[] psum = new double[ndim];
    double[] pmin = new double[ndim];
    double[] x = new double[ndim];
    p = new double[pp.length][pp[0].length];
    copyAssign(p,pp); // XXX p=pp, in NR "=" is overloading, see "nr3.h"
    //y = resize(y,mpts);
    y = new double[mpts];
    for (int i=0;i<mpts;i++) {
      for (int j=0;j<ndim;j++)
        x[j]=p[i][j];
      y[i]=func.funk(x);
    }
    nfunc=0;
    get_psum(p,psum);
    for (;;) {
      ilo=0;
    // First we must determine which point is the highest (worst), next-highest, and
    // lowest (best), by looping over the points in the simplex  
      // XXX ihi = y[0]>y[1] ? (inhi=1,0) : (inhi=0,1);
      ihi = 0;
      if(y[0]>y[1]) {
        inhi=1;
        ihi=0;
      }
      else {
        inhi=0;
        ihi=1;
      }
      
      for (int i=0;i<mpts;i++) {
        if (y[i] <= y[ilo]) ilo=i;
        if (y[i] > y[ihi]) {
          inhi=ihi;
          ihi=i;
        } else if (y[i] > y[inhi] && i != ihi) inhi=i;
      }
      double rtol=2.0*abs(y[ihi]-y[ilo])/(abs(y[ihi])+abs(y[ilo])+TINY);
      // Compute the fractional range from highest to lowest an   return if satisfactory
      if (rtol < ftol) {   // If returning, put best point and value in slot 0.
        swap(y,0,ilo);
        for (int i=0;i<ndim;i++) {
          // SWAP(p[0][i],p[ilo][i]);
          double dum = p[0][i]; p[0][i] = p[ilo][i]; p[ilo][i]=dum;
          pmin[i]=p[0][i];
        }
        fmin=y[0];
        return pmin;
      }
      if (nfunc >= NMAX) {
        throw new Error("NMAX exceeded");
      }
      nfunc += 2;
      // Begin a new iteration. First extrapolate by a factor -1 through the face of the 
      // simplex across from the high point, i.e. reflect the simplex from the high point.
      double ytry=amotry(p,y,psum,ihi,-1.0,func);
      if (ytry <= y[ilo])
          // Gives a result better than the best point, so try an additional extrapolation 
          // by a factor 2
        ytry=amotry(p,y,psum,ihi,2.0,func);
      else if (ytry >= y[inhi]) {
          // The reflected point is worse than the second-highest, so look for an intermediate 
          // lower point, i.e., do a one-dimensional contraction.
        double ysave=y[ihi];
        ytry=amotry(p,y,psum,ihi,0.5,func);
        if (ytry >= ysave) {    // Can't seem to get rid of that igh point.
          for (int i=0;i<mpts;i++) {  // Better contract around the lowest 
            if (i != ilo) {                    // (best) point
              for (int j=0;j<ndim;j++)
                p[i][j]=psum[j]=0.5*(p[i][j]+p[ilo][j]);
              y[i]=func.funk(psum);
            }
          }
          nfunc += ndim;    // Keep track of function ealuations
          get_psum(p,psum);    // Recompute psum
        }
      } else --nfunc;    // Correct the evaluation count.
    }    // Go back for the test of doneness and the next iteration
  }
  
  // Utility function
  public void get_psum(final double[][]p, final double[] psum) {
    for (int j=0;j<ndim;j++) {
      double sum=0.0;
      for (int i=0;i<mpts;i++)
        sum += p[i][j];
      psum[j]=sum;
    }
  }

  // Helper function: Extrapolates by a factor fac through the face of the simplex
  // across from the high point, tries it, and replaces the high point if the new point is better
  public double amotry(final double[][] p, final double[] y, final double[] psum,
    final int ihi, final double fac, final RealValueFun func) {
    double[] ptry = new double[ndim];
    double fac1=(1.0-fac)/ndim;
    double fac2=fac1-fac;
    for (int j=0;j<ndim;j++)
      ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
    double ytry=func.funk(ptry);    // Evaluate the function at the trial point.
    if (ytry < y[ihi]) {     // If it's better than the highest, then replace the 
      y[ihi]=ytry;           // highest.
      for (int j=0;j<ndim;j++) {
        psum[j] += ptry[j]-p[ihi][j];
        p[ihi][j]=ptry[j];
      }
    }
    return ytry;
  }
}
