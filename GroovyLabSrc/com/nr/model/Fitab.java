package com.nr.model;

import static com.nr.NRUtil.*;
import static java.lang.Math.*;
import com.nr.sf.Gamma;


/**
 * Fitting Data to a Straight Line
 * Copyright (C) Numerical Recipes Software 1986-2007
 * Java translation Copyright (C) Huang Wen Hui 2012
 *
 * @author hwh
 *
 */

/*
 Object for fitting a straight line y = a+bx to a set of points (x_i, y_i), with or without available
 errors σ_i. Call one of the two constructors to calculate the fit. The answers are then available 
 as the variables  a, b, siga, sigb, chi2 and either q or sigdat
 */
public class Fitab {
  int ndata;
  public double a, b, siga, sigb, chi2, q, sigdat;
  public double[] x, y, sig;

  /* 
   Constructor. Given a set of data points x[0..ndata-1], y[0..ndata-1] with individual standard 
   deviations sig[0..ndata-1], sets a, b and their respective probable uncertainties siga and sigb, the chi-square 
   chi2, and the goodness-of-fit probability q (that the fit would have χ^2 this large or larger)
   */
  public Fitab(final double[] xx, final double[] yy, final double[] ssig)  {
    ndata = xx.length;
    x = xx;
    y = yy;
    sig = ssig;
    chi2 = 0.;
    q = 1.0;
    sigdat = 0.;
    
    Gamma gam = new Gamma();
    int i;
    double ss=0.,sx=0.,sy=0.,st2=0.,t,wt,sxoss;
    b=0.0;   // accumulate sums ...
    for (i=0;i<ndata;i++) {
      wt=1.0/SQR(sig[i]);   // .. with weights
      ss += wt;
      sx += x[i]*wt;
      sy += y[i]*wt;
    }
    sxoss=sx/ss;
    for (i=0;i<ndata;i++) {
      t=(x[i]-sxoss)/sig[i];
      st2 += t*t;
      b += t*y[i]/sig[i];
    }
    b /= st2;    // Solve for α, b, σ_α, and σ_b
    a=(sy-sx*b)/ss;
    siga=sqrt((1.0+sx*sx/(ss*st2))/ss);
    sigb=sqrt(1.0/st2);     // Calculate χ^2
    for (i=0;i<ndata;i++) chi2 += SQR((y[i]-a-b*x[i])/sig[i]);
    if (ndata>2) q=gam.gammq(0.5*(ndata-2),0.5*chi2);
  }

  // Constructor. As above, but without known errors (sig is not used). The uncertainties siga
  // and sigb are estimated by assuming equal errors for all points, and that a straight line is 
  // a good fit. q is returned as 1.0, the normalization of chi2 is to unit standard deviation on all
  // points, and sigdat is set to the estimated error of each point.
  public Fitab(final double[] xx, final double[] yy) {
    ndata = xx.length;
    x = xx;
    y = yy;
    sig = xx;
    chi2 = 0.;
    q = 1.0;
    sigdat = 0.;
    int i;
    double ss,sx=0.,sy=0.,st2=0.,t,sxoss;
    b=0.0;    // Accumulate sums ...
    for (i=0;i<ndata;i++) {
      sx += x[i];     // .. without weights
      sy += y[i];
    }
    ss=ndata;
    sxoss=sx/ss;
    for (i=0;i<ndata;i++) {
      t=x[i]-sxoss;
      st2 += t*t;
      b += t*y[i];
    }
    b /= st2;       // Solve for α, b, σ_α, and σ_b
    a=(sy-sx*b)/ss;
    siga=sqrt((1.0+sx*sx/(ss*st2))/ss);
    sigb=sqrt(1.0/st2);    // Calclate χ^2
    for (i=0;i<ndata;i++) chi2 += SQR(y[i]-a-b*x[i]);
    if (ndata > 2) sigdat=sqrt(chi2/(ndata-2));    // For unweighted data evaluate typical sig using chi2, and adjust the standard deviations
    siga *= sigdat;
    sigb *= sigdat;
  }
}
