package com.nr.model;
import static com.nr.NRUtil.*;
import com.nr.RealMultiValueFun;
import com.nr.UniVarRealMultiValueFun;
import com.nr.la.SVD;

/**
 * general linear fit using SVD
 * 
 * Copyright (C) Numerical Recipes Software 1986-2007
 * Java translation Copyright (C) Huang Wen Hui 2012
 *
 * @author hwh
 *
 */

/*
 Object for general linear least-squares fitting using singular value decomposition. Call
constructor to bind data vectors and fitting functions. Then call fit, which sets
the output quantities a, covar, and chisq
 */
public class Fitsvd {
  int ndat, ma;
  final double tol;
  //double[] *x,&y,&sig;
  double[] x, y, sig;
  UniVarRealMultiValueFun funcs;
  public double[] a;  // Output values. a is the vector of fitted coefficients
  public double[][] covar;    // covar is its covariance matrix, and chisq is the value of χ^2 for the fit
  public double chisq;
  
  double[][] xmd;
  RealMultiValueFun funcsmd;
  
  // Constructor . Binds references to the data arrays xx, yy, and ssig, and to a user-supplied 
  // function funks(x) that returns a VecDoub containing ma basis functions evaluated at x = x.
  // If TOL is positive, it is the threshold (relative to the largest singular value) for discarding
  // small singular values. If it is <= 0, the default value in SVD is used
  public Fitsvd(final double[] xx, final double[] yy, final double[] ssig,
      final UniVarRealMultiValueFun funks) {
    this(xx,yy,ssig, funks, 1.e-12);
  }
  
  public Fitsvd(final double[] xx, final double[] yy, final double[] ssig,
      final UniVarRealMultiValueFun funks, final double TOL) {
    ndat = yy.length;
    x = xx;
    xmd = null;
    y = yy;
    sig = ssig;
    funcs = funks;
    tol = TOL;
  }

  // Solve by singular value decomposition the χ^2 minimization that fits for the coefficients 
  // a[0..ma-1] of a function that depends linearly on a, y = Σ_i a_i X funks_i(x). Set answer 
  // values for a[0..ma-1], chisq = χ^2, and the covariance matrix covar[0..ma-1][0..ma-1]
  public void fit() {
    int i,j,k;
    double tmp,thresh,sum;
    if (x!= null) ma = funcs.funk(x[0]).length;
    else ma = funcsmd.funk(row(xmd,0)).length;  // Discussed in 15.4.4.
    a = new double[ma];
    covar = new double[ma][ma];
    double[][] aa = new double[ndat][ma];
    double[] b = new double[ndat],afunc = new double[ma];
    for (i=0;i<ndat;i++) {   // Accumulate coefficients of the design matrix
      if (x!= null) afunc=funcs.funk(x[i]);
      else afunc=funcsmd.funk(row(xmd,i));   // Discussed in 15.4.4
      tmp=1.0/sig[i];
      for (j=0;j<ma;j++) aa[i][j]=afunc[j]*tmp;
      b[i]=y[i]*tmp;
    }
    SVD svd = new SVD(aa);   // Singular value decomposition
    thresh = (tol > 0. ? tol*svd.w[0] : -1.);
    svd.solve(b,a,thresh);    // Solve for the coefficients
    chisq=0.0;   // Evaluate chi-square
    for (i=0;i<ndat;i++) {
      sum=0.;
      for (j=0;j<ma;j++) sum += aa[i][j]*a[j];
      chisq += SQR(sum-b[i]);
    }
    for (i=0;i<ma;i++) {    // Sum contributions to covariance matrix (15.4.20)
      for (j=0;j<i+1;j++) {
        sum=0.0;
        for (k=0;k<ma;k++) if (svd.w[k] > svd.tsh)
          sum += svd.v[i][k]*svd.v[j][k]/SQR(svd.w[k]);
        covar[j][i]=covar[i][j]=sum;
      }
    }
  }


  public Fitsvd(final double[][] xx, final double[] yy, final double[] ssig,
      final RealMultiValueFun funks) {
    this(xx,yy,ssig,funks,1.e-12);
  }
  
  // Constructor for multidimensional fits. Exactly the same as the previous constructor,
  // except that xx is now a matrix whose rows are the multidimensional data points and funks is
  // now a function of a multidimensional data point (as a VecDoub)
  public Fitsvd(final double[][] xx, final double[] yy, final double[] ssig,
      final RealMultiValueFun funks, final double TOL) {
    ndat = yy.length;
    x = null;
    xmd = xx;
    y = yy;
    sig = ssig;
    funcsmd = funks;
    tol = TOL;
  }

  public double[] row(final double[][] a, final int i) {
    int j,n=a[0].length;
    double[] ans = new double[n];
    for (j=0;j<n;j++) ans[j] = a[i][j];
    return ans;
  }
}
