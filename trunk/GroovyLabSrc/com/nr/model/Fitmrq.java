package com.nr.model;

import static com.nr.NRUtil.*;
import static com.nr.la.GaussJordan.*;
import static java.lang.Math.*;

import org.netlib.util.doubleW;


/**
 * Levenberg-Marquardt nonlinear fitting
 * Copyright (C) Numerical Recipes Software 1986-2007
 * Java translation Copyright (C) Huang Wen Hui 2012
 *
 * @author hwh
 *
 */
/*
 Object for nonlinear least-squares fitting by the Levenberg-Marquardt method, also 
 including the ability to hold specified parameters at fixed, specified values. Call constructor
 to bind data vectors and fitting functions and to input an initial parameter guess. Then call
 any combination of hold, free, and fit as often as desired. fit sets the output quantities
 a, covar, alpha, and chisq
 */
public class Fitmrq {
  static final int NDONE=4, ITMAX=1000;
  int ndat, ma, mfit;
  double[] x,y,sig;
  final double tol;
  MultiFuncd funcs;
  public boolean[] ia;
  public double[] a;
  public double[][] covar;
  public double[][] alpha;
  public double chisq;

  /*
   Constructor. Binds references to the data arrays xx, yy, and ssig, and to a user-supplied function
   funks that calculates the nonlinear fitting function and its derivatives. Also inputs the initial 
   parameters guess aa (which is copied, not modified) and an optional convergence 
   tolerance TOL. Initializes all parameters as free (not held)
   */
  public Fitmrq(final double[] xx, final double[] yy, final double[] ssig, final double[] aa,
      final MultiFuncd funks) {
    this(xx, yy, ssig, aa, funks, 1.e-3);
  }
  
  public Fitmrq(final double[] xx, final double[] yy, final double[] ssig, final double[] aa,
      final MultiFuncd funks, final double TOL) {
    ndat = xx.length;
    ma = aa.length;
    x = xx;
    y = yy;
    sig = ssig;
    tol = TOL;
    funcs = funks;
    ia = new boolean[ma];
    alpha = new double[ma][ma];
    a = buildVector(aa);
    covar = new double[ma][ma];
    for (int i=0;i<ma;i++) ia[i] = true;
  }

  /*
   Optional functions for holding a parameter, identified by a value i in the range 0, .., ma-1,
   fized at the value val, or for freeing a parameter that was previously held fixed.
   hold and free may be called for any number of parameters before calling fit to calculate best-fit
   values for the remaining (not held) parameters, and the process may be repeated multiple times.
   */
  public void hold(final int i, final double val) {ia[i]=false; a[i]=val;}
  
  public void free(final int i) {ia[i]=true;}

  /*
   Iterate to reduce the χ^2 of a fit between a set of data points x[0..ndat-1],
   y[0..ndat-1] with individual standard deviations sig[0..ndat-1], and a nonlinear
   function that depends on ma coefficients a[0..ma-1]. When χ^2 is no longer 
   decreasing, set best-fit values for the parameters a[0..ma-1], and chisq = χ^2,
   covar[0..ma-1][0..ma-1], and alpha[0..ma-1][0..ma-1] (Parameters held fixed will return
   zero covariances).
   */
  public void fit() {
    int j,k,l,iter,done=0;
    double alamda=.001,ochisq;
    double[] atry = new double[ma];
    double[] beta = new double[ma];
    double[] da = new double[ma];
    ;
    mfit=0;
    for (j=0;j<ma;j++) if (ia[j]) mfit++;
    double[][] oneda = new double[mfit][1], temp = new double[mfit][mfit];
    mrqcof(a,alpha,beta);    // Initialization
    for (j=0;j<ma;j++) atry[j]=a[j];
    ochisq=chisq;
    for (iter=0;iter<ITMAX;iter++) {
      if (done==NDONE) alamda=0.;   // Last pass. Use zero alamda
      for (j=0;j<mfit;j++) {   // After linearized fitting matrix, by agmenting diagonal elements
        for (k=0;k<mfit;k++) covar[j][k]=alpha[j][k];
        covar[j][j]=alpha[j][j]*(1.0+alamda);
        for (k=0;k<mfit;k++) temp[j][k]=covar[j][k];
        oneda[j][0]=beta[j];
      }
      gaussj(temp,oneda);    // Matrix solution
      for (j=0;j<mfit;j++) {
        for (k=0;k<mfit;k++) covar[j][k]=temp[j][k];
        da[j]=oneda[j][0];
      }
      if (done==NDONE) {   // Converged. Clean up and return.
        covsrt(covar);
        covsrt(alpha);
        return;
      }
      for (j=0,l=0;l<ma;l++)    // Did the trial succeed?
        if (ia[l]) atry[l]=a[l]+da[j++];
      mrqcof(atry,covar,da);
      if (abs(chisq-ochisq) < max(tol,tol*chisq)) done++;
      if (chisq < ochisq) {    // Success, accept the new solution.
        alamda *= 0.1;
        ochisq=chisq;
        for (j=0;j<mfit;j++) {
          for (k=0;k<mfit;k++) alpha[j][k]=covar[j][k];
            beta[j]=da[j];
        }
        for (l=0;l<ma;l++) a[l]=atry[l];
      } else {   // Failure, increase alamda.
        alamda *= 10.0;
        chisq=ochisq;
      }
    }
    throw new IllegalArgumentException("Fitmrq too many iterations");
  }


  // Used by fit to evaluate the linearized fitting matrix alpha, and vector beta as in (15.5.8) and
  // to calculate χ^2
  public void mrqcof(final double[] a, final double[][] alpha, final double[] beta) {
    int i,j,k,l,m;
    double ymod,wt,sig2i,dy;
    doubleW ymodW= new doubleW(0);
    double[] dyda = new double[ma];
    for (j=0;j<mfit;j++) {   // Initialize (symmetric) alpha, beta
      for (k=0;k<=j;k++) alpha[j][k]=0.0;
      beta[j]=0.;
    }
    chisq=0.;
    for (i=0;i<ndat;i++) {   // Summation loop over all data
      funcs.funk(x[i],a,ymodW,dyda); ymod = ymodW.val;
      sig2i=1.0/(sig[i]*sig[i]);
      dy=y[i]-ymod;
      for (j=0,l=0;l<ma;l++) {
        if (ia[l]) {
          wt=dyda[l]*sig2i;
          for (k=0,m=0;m<l+1;m++)
            if (ia[m]) alpha[j][k++] += wt*dyda[m];
          beta[j++] += dy*wt;
        }
      }
      chisq += dy*dy*sig2i;   // And find χ^2
    }
    for (j=1;j<mfit;j++)  //  Fill in the symmetric side
      for (k=0;k<j;k++) alpha[k][j]=alpha[j][k];
  }

  // Expand in storage the covariance matrix covar, so as to take into account parameters that are being 
  // held fixed (For the later, return zero covariances.)
  public void covsrt(final double[][] covar) {
    int i,j,k;
    for (i=mfit;i<ma;i++)
      for (j=0;j<i+1;j++) covar[i][j]=covar[j][i]=0.0;
    k=mfit-1;
    for (j=ma-1;j>=0;j--) {
      if (ia[j]) {
        for (i=0;i<ma;i++) {
          //SWAP(covar[i][k],covar[i][j]);
          double swap = covar[i][k]; covar[i][k] = covar[i][j]; covar[i][j] = swap;
        }
        for (i=0;i<ma;i++) {
          //SWAP(covar[k][i],covar[j][i]);
          double swap = covar[k][i]; covar[k][i] = covar[j][i]; covar[j][i] = swap;
        }
        k--;
      }
    }
  }

}
