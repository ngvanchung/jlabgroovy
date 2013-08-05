package com.nr.ci;
import static java.lang.Math.*;
import static com.nr.NRUtil.*;
import com.nr.ran.Ran;
import com.nr.sort.Indexx;

/**
 * Support Vector Machines
 * Copyright (C) Numerical Recipes Software 1986-2007
 * Java translation Copyright (C) Huang Wen Hui 2012
 *
 * @author hwh
 *
 */

// class for solving SVM problems by the SOR method
public class Svm {
  private Svmgenkernel gker;  // Reference bound to user's kernel (and data)
  private int m, fnz, fub, niter;
  private double[] alph, alphold;   // Vectors of a's before and after a step
  private Ran ran;   // Random number generator
  private boolean alphinit;
  private double dalph;   // Change in norm of the a's in one step
  
  // constructor binds the user's kernel and allocates storage
  public Svm(final Svmgenkernel inker){
    gker = inker;
    m = gker.y.length;
    alph = new double[m];
    alphold= new double[m];
    ran = new Ran(21);
    alphinit = false;
  }
  
  // Perform one group of relaxation steps: a single step over all the a's, and multiple steps over only the interior a's
  public double relax(final double lambda, final double om) {
    int iter,j,jj,k,kk;
    double sum;   // index when a's are sorted by value
    double[] pinsum = new double[m];  // stored sums over noninterior variables
    if (alphinit == false) {   // start all a's at 0
      for (j=0; j<m; j++) alph[j] = 0.;
      alphinit = true;
    }
    alphold = alph;   // save old a's
    // here begins the relaxation pass over all the a's
    Indexx x = new Indexx(alph);  // sort a's, then find first nonzero one
    for (fnz=0; fnz<m; fnz++) if (alph[x.indx[fnz]] != 0.) break; 
    for (j=fnz; j<m-2; j++) {  // randomly permute all the nonzero a's
      k = j + (ran.int32p() % (m-j));
      swap(x.indx,j,k);
    }
    for (jj=0; jj<m; jj++) {   // Main loop over a's
      j = x.indx[jj];
      sum = 0.;
      for (kk=fnz; kk<m; kk++) {  // Sums start with first nonzero
        k = x.indx[kk];
        sum += (gker.ker[j][k] + 1.)*gker.y[k]*alph[k];
      }
      alph[j] = alph[j] - (om/(gker.ker[j][j]+1.))*(gker.y[j]*sum-1.);
      alph[j] = max(0.,min(lambda,alph[j]));   // Projection operator
      if (jj < fnz && alph[j]!=0) {
        --fnz;
        //SWAP(x.indx[--fnz],x.indx[jj]);
        swap(x.indx, fnz, jj);
      }  // (Above) MAke an \alpha active if it becomes nonzero
    }
    
    // Here begins the relaxation passes over the interior \alpha's
    Indexx y = new Indexx(alph);   // Sort. Identify interior \alpha's
    for (fnz=0; fnz<m; fnz++) if (alph[y.indx[fnz]] != 0.) break; 
    for (fub=fnz; fub<m; fub++) if (alph[y.indx[fub]] == lambda) break;
    for (j=fnz; j<fub-2; j++) {  // Permute
      k = j + (ran.int32p() % (fub-j));
      swap(y.indx,j,k);
    }
    
    for (jj=fnz; jj<fub; jj++) {  // Compute sums over pinned \alpha's just once
      j = y.indx[jj];
      sum = 0.;
      for (kk=fub; kk<m; kk++) {
        k = y.indx[kk];
        sum += (gker.ker[j][k] + 1.)*gker.y[k]*alph[k];
      }
      pinsum[jj] = sum;
    }
    niter = max((int)(0.5*(m+1.0)*(m-fnz+1.0)/(SQR(fub-fnz+1.0))),1);
    
    // Calculate a numer of iterations that will take about half as long as the full pass just completed
    for (iter=0; iter<niter; iter++) {  // Main loop over \alpha's
      for (jj=fnz; jj<fub; jj++) {
        j = y.indx[jj];
        sum = pinsum[jj];
        for (kk=fnz; kk<fub; kk++) {
          k = y.indx[kk];
          sum += (gker.ker[j][k] + 1.)*gker.y[k]*alph[k];
        }
        alph[j] = alph[j] - (om/(gker.ker[j][j]+1.))*(gker.y[j]*sum-1.);
        alph[j] = max(0.,min(lambda,alph[j]));
      }   
    }
    dalph = 0.;  // Return change in norm of \alpha vector
    for (j=0;j<m;j++) dalph += SQR(alph[j]-alphold[j]);
    return sqrt(dalph);
  }
  
  
  // Call only after convergence via repeated calls to relax. Returns the decision rule f(x) for data point k
  public double predict(final int k) {
    double sum = 0.;
    for (int j=0; j<m; j++) sum += alph[j]*gker.y[j]*(gker.ker[j][k]+1.0);
    return sum;
  }
  
  // Call only after convergence via repeated calls to relax. Returns the decision rule f(x) for an arbitrary feature vector
  public double predict(final double[] x) {
    double sum = 0.;
    for (int j=0; j<m; j++) sum += alph[j]*gker.y[j]*(gker.kernel(j,x)+1.0);
    return sum;
  }
}
