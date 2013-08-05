package com.nr.sf;
import static java.lang.Math.*;
import static com.nr.NRUtil.*;

// Logistic distribution
public class Logisticdist {
  double mu, sig;
  // Logistic(0, 1)
  public Logisticdist(){
    this(0.,1.);
  }
  
  public Logisticdist(final double mmu, final double ssig) {
    mu = mmu;
    sig=ssig;
    if (sig <= 0.) throw new IllegalArgumentException("bad sig in Logisticdist");
  }
  
  // probability density function
  public double p(final double x) {
    double e = exp(-abs(1.81379936423421785*(x-mu)/sig));
    return 1.81379936423421785*e/(sig*SQR(1.+e));
  }
  
  // cumulative distribution function
  public double cdf(final double x) {
    double e = exp(-abs(1.81379936423421785*(x-mu)/sig));
    // because we used abs to control overflow, we now have two cases
    if (x >= mu) return 1./(1.+e);
    else return e/(1.+e);
  }
  
  // inverse cumulative distribution function
  public double invcdf(final double p) {
    if (p <= 0. || p >= 1.) throw new IllegalArgumentException("bad p in Logisticdist");
    return mu + 0.551328895421792049*sig*log(p/(1.-p));
  }
}
