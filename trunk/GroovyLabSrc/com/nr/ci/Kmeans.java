package com.nr.ci;
import static com.nr.NRUtil.*;

/**
 * k-means classification
 * Copyright (C) Numerical Recipes Software 1986-2007
 * Java translation Copyright (C) Huang Wen Hui 2012
 *
 * @author hwh
 *
 */

/*
 Solve for a k-means clustering model from a set of data points and initial guesses of the means.
 Output is a set of means and an assignment of each data pointto one component.
 */
public class Kmeans {
  public int nn, mm, kk, nchg;
  public double[][] data, means;
  public int[] assign, count;
  
  /*
   Arguments are the data points (as rows in a matrix), and initial guesses for the means (also as rows in a matrix)
   */
  public Kmeans(final double[][] ddata, final double[][] mmeans) {
    nn = ddata.length;
    mm = ddata[0].length;
    kk = mmeans.length;
    data = buildMatrix(ddata);
    means = buildMatrix(mmeans);
    assign = new int[nn];
    count = new int[kk];
    
    /*
     Perform one initial -step and M-step. User is responsible for calling additional steps until 
     convergence is obtained.
     */
    estep();
    mstep();
  }
  
  public int estep() {
    int k,m,n,kmin=0;
    double dmin,d;
    nchg = 0;
    for (k=0;k<kk;k++) count[k] = 0;
    for (n=0;n<nn;n++) {
      dmin = 9.99e99;
      for (k=0;k<kk;k++) {
        for (d=0.,m=0; m<mm; m++) d += SQR(data[n][m]-means[k][m]);
        if (d < dmin) {dmin = d; kmin = k;}
      }
      if (kmin != assign[n]) nchg++;
      assign[n] = kmin;
      count[kmin]++;
    }
    return nchg;
  }
  
  public void mstep() {
    int n,k,m;
    for (k=0;k<kk;k++) for (m=0;m<mm;m++) means[k][m] = 0.;
    for (n=0;n<nn;n++) for (m=0;m<mm;m++) means[assign[n]][m] += data[n][m];
    for (k=0;k<kk;k++) {
      if (count[k] > 0) for (m=0;m<mm;m++) means[k][m] /= count[k];
    }
  }
}
