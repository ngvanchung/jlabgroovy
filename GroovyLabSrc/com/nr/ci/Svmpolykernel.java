package com.nr.ci;

import static java.lang.Math.*;

// kernel structure for the polynomial kernel
public class Svmpolykernel extends Svmgenkernel {
  int n;

  double a, b, d;

  // constructor is called with the mXn data matrix, the vector of y_i's length m, and the 
  // constants aa, bb, and dd
  public Svmpolykernel(final double[][] ddata, final double[] yy,
      final double aa, final double bb, final double dd) {
    super(yy, ddata);
    n = data[0].length;
    a = aa;
    b = bb;
    d = dd;
    fill();
  }

  public double kernel(final double[] xi, final double[] xj) {
    double dott = 0.;
    for (int k = 0; k < n; k++)
      dott += xi[k] * xj[k];
    return pow(a * dott + b, d);
  }
}
