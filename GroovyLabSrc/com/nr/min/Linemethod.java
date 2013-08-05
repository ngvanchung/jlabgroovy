package com.nr.min;

import com.nr.RealValueFun;
// Base class for line-minimization algorithms. Provides the line-minimization routine linmin
public class Linemethod {
  public double[] p;
  public double[] xi;
  RealValueFun func;
  int n;
  
  // Constructor argument is the user-supplied function to be minimized
  public Linemethod(final RealValueFun funcc) {
    func = funcc;
  }
  
  
  /*
   Line-minimization routine. Given an n-dimensional point p[0..n-1] and an n-dimensional
   direction xi[0..n-1], moves and resets p to where the function func(p) takes on a 
   minimum along the direction xi from p, and replaces xi by the actual vector displacement\
   that p was moved. Also returns the value of func at the returned location p. This is 
   actually all accomplished by calling the routine bracket and minimize of Brent.
   */
  public double linmin() {
    double ax,xx,xmin;
    n=p.length;
    F1dim f1dim = new F1dim(p,xi,func);
    ax=0.0;    // Initial guess for brackets
    xx=1.0;
    Brent brent = new Brent();
    brent.bracket(ax,xx,f1dim);
    xmin=brent.minimize(f1dim);
    for (int j=0;j<n;j++) {  // Construct the vector results and return
      xi[j] *= xmin;
      p[j] += xi[j];
    }
    return brent.fmin;
  }
}
