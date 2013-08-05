package com.nr.model;

import static com.nr.NRUtil.*;
import static java.lang.Math.*;
import static com.nr.sort.Sorter.*;
/**
 * fit a line minimizing absolute deviation
 * 
 * Copyright (C) Numerical Recipes Software 1986-2007
 * Java translation Copyright (C) Huang Wen Hui 2012
 *
 * @author hwh
 *
 */

/*
 Object for fitting a straight line y = a+b*x to a set of points (x_i, y_i), by the criterion of least 
 absolute deviations. Call the constructor to calculate the fit. The answers are then available as the variables a, b, and 
 abdev (the mean absolute deviation of the points from the line).
 */
public class Fitmed {
  int ndata;
  public double a, b, abdev;
  double[] x, y;

  // Constructor. Given a set of data points xx[0..ndata-1], yy[0..ndata-1], sets a, b, and abdev
  public Fitmed(final double[] xx, final double[] yy) {
    ndata = xx.length;
    x = xx;
    y = yy;
    
    int j;
    double b1,b2,del,f,f1,f2,sigb,temp;
    double sx=0.0,sy=0.0,sxy=0.0,sxx=0.0,chisq=0.0;
    for (j=0;j<ndata;j++) {  // As a first guess for a and b, we will find the least-squares fitting line
      sx += x[j];
      sy += y[j];
      sxy += x[j]*y[j];
      sxx += SQR(x[j]);
    }
    del=ndata*sxx-sx*sx;
    a=(sxx*sy-sx*sxy)/del;   // Least-squares solutions.
    b=(ndata*sxy-sx*sy)/del;
    for (j=0;j<ndata;j++) {
      temp=y[j]-(a+b*x[j]);
      chisq += (temp*temp);
    }
    sigb=sqrt(chisq/del);   // The standard deviation will give some idea of how big an iteration step to take.
    b1=b;
    f1=rofunc(b1);
    if (sigb > 0.0) {
      b2=b+SIGN(3.0*sigb,f1);   // Guess bracket as 3-σ away, in the downhill direction known from f1
      f2=rofunc(b2);
      if (b2 == b1) {
        abdev /= ndata;
        return;
      }
      while (f1*f2 > 0.0) {   // Bracketing
        b=b2+1.6*(b2-b1);
        b1=b2;
        f1=f2;
        b2=b;
        f2=rofunc(b2);
      }
      sigb=0.01*sigb;
      while (abs(b2-b1) > sigb) {
        b=b1+0.5*(b2-b1);    // Bisection
        if (b == b1 || b == b2) break;
        f=rofunc(b);
        if (f*f1 >= 0.0) {
          f1=f;
          b1=b;
        } else {
          f2=f;
          b2=b;
        }
      }
    }
    abdev /= ndata;
  }

  // Evaluates the right-hand side of equation (15.7.16) for a given value of b
  public double rofunc(final double b) {
    final double EPS=DBL_EPSILON;
    int j;
    double d,sum=0.0;
    double[] arr = new double[ndata];
    for (j=0;j<ndata;j++) arr[j]=y[j]-b*x[j];
    if ((ndata & 1) == 1) {
      a=select((ndata-1)>>1,arr);
    } else {
      j=ndata >> 1;
      a=0.5*(select(j-1,arr)+select(j,arr));
    }
    abdev=0.0;
    for (j=0;j<ndata;j++) {
      d=y[j]-(b*x[j]+a);
      abdev += abs(d);
      if (y[j] != 0.0) d /= abs(y[j]);
      if (abs(d) > EPS) sum += (d >= 0.0 ? x[j] : -x[j]);
    }
    return sum;
  }
}
