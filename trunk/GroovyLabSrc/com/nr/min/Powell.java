package com.nr.min;

import static com.nr.NRUtil.*;
import static java.lang.Math.*;
import com.nr.RealValueFun;

/* Multidimensional minimization by Powell's method */
// Constructor arguments are func the function to be minimized,
// and an optional argument ftoll, the fractional tolerance in the function value such
// that failure to decrease by more than this amount on one iteration signals doneness

public class Powell extends Linemethod{
  public int iter;
  public double fret;
  
  /*
  using Linemethod<T>::func;
  using Linemethd<T>::linmin;
  using Linemethod<T>::p;
  using Linemethod<T>::xi;
  */
  
  
  final double ftol;
  
  public Powell(final RealValueFun func){
    this(func, 3.0e-8);
  }
  
  public Powell(final RealValueFun func, final double ftoll) {
    super(func);
    ftol = ftoll;
  }
  
// minimization of a function of n variables. Input consists of an initial starting point pp[0..n-1]
// Returned is the best point found, at which point fret
//  is the minimum function value and iter is the number of iterations taken 
  public double[] minimize(final double[] pp) {
    int n=pp.length;
    double[][] ximat = new double[n][n];
    for (int i=0;i<n;i++) ximat[i][i]=1.0;
    return minimize(pp,ximat);
  }
  
  // alternative interface: Input consists of the initial starting point pp[0..n-1] and an initial
  // matrix ximat[0..n-1][0..n-1], whose columns contain the initial set of directions.
  // On output ximat is the then-current direction set
  public double[] minimize(final double[] pp, final double[][] ximat) {
    final int ITMAX=200;
    final double TINY=1.0e-25;
    double fptt;
    int n=pp.length;
    p = buildVector(pp);
    double[] pt= new double[n], ptt = new double[n];
    // xi = resize(xi,n);
    xi = new double[n];
    fret=func.funk(p);
    for (int j=0;j<n;j++) pt[j]=p[j];
    for (iter=0;;++iter) {
      double fp=fret;
      int ibig=0;
      double del=0.0;
      for (int i=0;i<n;i++) {
        for (int j=0;j<n;j++) xi[j]=ximat[j][i];
        fptt=fret;
        fret=linmin();
        if (fptt-fret > del) {
          del=fptt-fret;
          ibig=i+1;
        }
      }
      if (2.0*(fp-fret) <= ftol*(abs(fp)+abs(fret))+TINY) {
        return p;
      }
      if (iter == ITMAX) throw new IllegalArgumentException("powell exceeding maximum iterations.");
      for (int j=0;j<n;j++) {
        ptt[j]=2.0*p[j]-pt[j];
        xi[j]=p[j]-pt[j];
        pt[j]=p[j];
      }
      fptt=func.funk(ptt);
      if (fptt < fp) {
        double t=2.0*(fp-2.0*fret+fptt)*SQR(fp-fret-del)-del*SQR(fp-fptt);
        if (t < 0.0) {
          fret=linmin();
          for (int j=0;j<n;j++) {
            ximat[j][ibig-1]=ximat[j][n-1];
            ximat[j][n-1]=xi[j];
          }
        }
      }
    }
  }
}
