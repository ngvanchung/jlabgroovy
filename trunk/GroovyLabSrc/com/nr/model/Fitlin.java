package com.nr.model;
import static com.nr.NRUtil.*;
import static com.nr.la.GaussJordan.*;
import com.nr.UniVarRealMultiValueFun;

/**
 * General linear fit
 * Copyright (C) Numerical Recipes Software 1986-2007
 * Java translation Copyright (C) Huang Wen Hui 2012
 *
 * @author hwh
 *
 */
/* 
 Object  for general linear least-squares fitting by solving the normal equations, also including
 the ability to hold specified parameters at fixed, specified values. Call constructor to bind data
 vectors and fitting functions. Then call any combination of hold, free, and fit as often
 as desired. fit sets the output quantities a, covar, and chisq
 */

public class Fitlin {
  int ndat, ma;
  public double[] x,y,sig;
  UniVarRealMultiValueFun funcs;
  public boolean[] ia;

  public double[] a;    // Output values. a is the vector of fitted coefficients
  public double[][] covar;   // covar is its covariance matrix, and chisq is
  public double chisq;   //  the value of χ^2 for the fit

  // Constructor. Binds references to the data arrays xx, yy, and ig, and to a user-supplied
  // function funcs(x) that returns a VecDoub containing ma basis functions evaluated at x = x
  // Initializes all parameters as free (not held)
  public Fitlin(final double[] xx, final double[] yy, final double[] ssig, final UniVarRealMultiValueFun funcs) {
    ndat = xx.length;
    x = xx;
    y = yy;
    sig = ssig;
    this.funcs = funcs;
    ma = funcs.funk(x[0]).length;
    a = new double[ma];
    covar = new double[ma][ma];
    ia = new boolean[ma];
    
    for (int i=0;i<ma;i++) ia[i] = true;    
  }

  // Optional functions for holding a parameter, identified by a value i in the range 0, ..., ma-1, fixed at the value
  // val, or for freeing a parameter that was previously held fixed. hold and free may be called for any number of
  // parameters before calling fit to calculate best-fit values for the remaining (not held|) parameters, and the process 
  // may be repeated multiple times. Alternatively, you can set the boolean vector ia directly, before calling fit
  public void hold(final int i, final double val) {ia[i]=false; a[i]=val;}
  public void free(final int i) {ia[i]=true;}
  
  /* 
   Solve the normal equations for χ^2 minimization to fit for some r all of the coefficients a[0..ma-1]
   of a function that depends linearly on a, y = Σ_i  a_i X funks_i(x). Set answer
   values for a[0..ma-1], χ^2 = chisq, and the covariance matrix covar[0..ma-1][0..ma-1]. 
   (Parameters held fixed by calls to hold will return zero covariances.)
   
   */
  public void fit() {
    int i,j,k,l,m,mfit=0;
    double ym,wt,sum,sig2i;
    double[] afunc = new double[ma];
    for (j=0;j<ma;j++) if (ia[j]) mfit++;
    if (mfit == 0) throw new IllegalArgumentException("lfit: no parameters to be fitted");
    double[][] temp = new double[mfit][mfit], beta = new double[mfit][1];
    for (i=0;i<ndat;i++) {    // Loop over data to accumulate coefficients of the normal equations
      afunc = funcs.funk(x[i]);
      ym=y[i];
      if (mfit < ma) {    // Subtract all dependences on known pieces of the fitting function
        for (j=0;j<ma;j++)
          if (!ia[j]) ym -= a[j]*afunc[j];
      }
      sig2i=1.0/SQR(sig[i]);
      for (j=0,l=0;l<ma;l++) {   // Set up matrix and r, h, s for matrix inversion
        if (ia[l]) {
          wt=afunc[l]*sig2i;
          for (k=0,m=0;m<=l;m++)
            if (ia[m]) temp[j][k++] += wt*afunc[m];
          beta[j++][0] += ym*wt;
        }
      }
    }
    for (j=1;j<mfit;j++) for (k=0;k<j;k++) temp[k][j]=temp[j][k];
    gaussj(temp,beta);    // Matrix solution
    for (j=0,l=0;l<ma;l++) if (ia[l]) a[l]=beta[j++][0];
     // Spread the solution to appropriate positions in a, and evaluate χ^2 of the fit
    chisq=0.0;
    for (i=0;i<ndat;i++) {
      afunc = funcs.funk(x[i]);
      sum=0.0;
      for (j=0;j<ma;j++) sum += a[j]*afunc[j];
      chisq += SQR((y[i]-sum)/sig[i]);
    }
    for (j=0;j<mfit;j++) for (k=0;k<mfit;k++) covar[j][k]=temp[j][k];
    for (i=mfit;i<ma;i++)  // Rearrange covariance matrix into correct order
      for (j=0;j<i+1;j++) covar[i][j]=covar[j][i]=0.0;
    k=mfit-1;
    for (j=ma-1;j>=0;j--) {
      if (ia[j]) {
        for (i=0;i<ma;i++){
          //SWAP(covar[i][k],covar[i][j]);
          double swap = covar[i][k]; covar[i][k] = covar[i][j]; covar[i][j] = swap;
        }
        for (i=0;i<ma;i++){
          // SWAP(covar[k][i],covar[j][i]);
          double swap = covar[k][i]; covar[k][i] = covar[j][i]; covar[j][i] =swap;
        }
        k--;
      }
    }
  }
}
