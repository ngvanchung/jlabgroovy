package com.nr.sf;

import static java.lang.Math.*;
import static com.nr.sf.Gamma.*;
import static com.nr.NRUtil.*;

// Student-t distribution derived from the beta function Beta
public class Studenttdist extends Beta {
  double nu, mu, sig, np, fac;
  
  // Initialize as Student(ν, 0, 1)
  public Studenttdist(double nnu) {
    this(nnu, 0., 1.);
  }
  
  // Initialize with ν, μ and σ. 
  public Studenttdist(double nnu, double mmu, double ssig ){
    nu = nnu;
    mu = mmu;
    sig = ssig;
    
    if (sig <= 0. || nu <= 0.) throw new IllegalArgumentException("bad sig,nu in Studentdist");
    np = 0.5*(nu + 1.);
    fac = gammln(np)-gammln(0.5*nu);
  }
  
  // probability density function
  public double p(double t) {
    return exp(-np*log(1.+SQR((t-mu)/sig)/nu)+fac)
      /(sqrt(3.14159265358979324*nu)*sig);
  }
  
  // cumulative distribution function
  public double cdf(double t) {
    double p = 0.5*betai(0.5*nu, 0.5, nu/(nu+SQR((t-mu)/sig)));
    if (t >= mu) return 1. - p;
    else return p;
  }
  
  // inverse cumulative distribution function
  public double invcdf(double p) {
    if (p <= 0. || p >= 1.) throw new IllegalArgumentException("bad p in Studentdist");
    double x = invbetai(2.*min(p,1.-p), 0.5*nu, 0.5);
    x = sig*sqrt(nu*(1.-x)/x);
    return (p >= 0.5? mu+x : mu-x);
  }
  
  // the two-sided cdf A(t|ν)
  public double aa(double t) {
    if (t < 0.) throw new IllegalArgumentException("bad t in Studentdist");
    return 1.-betai(0.5*nu, 0.5, nu/(nu+SQR(t)));
  }
  
  // the inverse, namely t such that p = A(t|ν)
  public double invaa(double p) {
    if (p < 0. || p >= 1.) throw new IllegalArgumentException("bad p in Studentdist");
    double x = invbetai(1.-p, 0.5*nu, 0.5);
    return sqrt(nu*(1.-x)/x);
  }
}
