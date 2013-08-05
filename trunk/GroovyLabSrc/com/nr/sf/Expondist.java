package com.nr.sf;
import static java.lang.Math.*;

// Exponential distribution
public class Expondist {
  double bet;
  
 // initialize with β
  public Expondist(final double bbet){
    bet = bbet;
    if (bet <= 0.) throw new IllegalArgumentException("bad bet in Expondist");  
  }
  
  // probability density function
  public double p(final double x) {
    if (x < 0.) throw new IllegalArgumentException("bad x in Expondist");
    return bet*exp(-bet*x);
  }
  
  // cumulative distribution function
  public double cdf(final double x) {
    if (x < 0.) throw new IllegalArgumentException("bad x in Expondist");
    return 1.-exp(-bet*x);
  }
  
  // inverse cumulative distribution function
  public double invcdf(final double p) {
    if (p < 0. || p >= 1.) throw new IllegalArgumentException("bad p in Expondist");
    return -log(1.-p)/bet;
  }
}
