package com.nr.min;

import com.nr.RealValueFun;
import com.nr.UniVarRealValueFun;

// Must accompany linmin in Linemethod
public class F1dim implements UniVarRealValueFun{
  final double[] p;
  final double[] xi;
  int n;
  RealValueFun func;
  double[] xt;
  
  // Constructor takes as inputs an n-dimensional point p[0..n-1] and an
  // n-dimensional direction xi[0..n-1] from linmin, as well as the function
  // that takes a vector argument
  public F1dim(final double[] pp, final double[] xii, final RealValueFun funcc) {
    p = pp;
    xi = xii;
    n = pp.length;
    func = funcc;
    xt = new double[n];
  }
  
  public double funk(final double x){
    return get(x);
  }
  
  public double get(final double x) {
    for (int j=0;j<n;j++)
      xt[j]=p[j]+x*xi[j];
    return func.funk(xt);
  }
}
