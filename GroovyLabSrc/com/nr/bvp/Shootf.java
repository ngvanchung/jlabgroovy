package com.nr.bvp;

import static com.nr.NRUtil.*;
import com.nr.RealMultiValueFun;
import com.nr.ode.DerivativeInf;
import com.nr.ode.Odeint;
import com.nr.ode.Output;
import com.nr.ode.StepperDopr853;

/* 
 For use with newt to solve a two-point boundary value problem by shooting to a fitting point.
 */
public abstract class Shootf implements RealMultiValueFun{
  int nvar,n2;
  double x1,x2,xf;
  double atol,rtol;
  double h1,hmin;
  double[] y,f1,f2;
  DerivativeInf d;
 
 /*
   Routine for use with newt to solve a two-point boundary value problem for nvar coupled ODEs
   by shooting from x1 and x2 to a fitting point xf. Initial values for the nvar ODEs at x1 are generated from the n2 
   coefficients v1 and the user-supplied routine load1. Likewise, those at x2 are from the n1 = nvar-n2 coeffcients v2,
   using load2. The coefficients v1 and v2 should be stored in a single array v[0..nvar-1] in the main program
   with v1 in v[0..n2-1] and v2 in v[n2..nvar-1]
   
   */
  public Shootf(final int nvarr, final int nn2, final double xx1, final double xx2, final double xxf, final DerivativeInf dd) {
    nvar=nvarr;
    n2=nn2;
    x1=xx1;
    x2=xx2;
    xf=xxf;
    d=dd;
    atol=1.0e-14;
    rtol = atol;
    hmin=0.0;
    y = new double[nvar];
    f1 = new double[nvar];
    f2 = new double[nvar];
  }
  
  /*
   This routine is used by newt. It integrates the ODEs to xf using an eighth-order Runge-Kutta method with absolute
   and relative tolerances atol and rtol, initial stepsize h1, and minimum stepsize hmin. At xf it calls the user-supplied 
   routine score to evaluate the nvar functions f1 and f2 that ougt to match at xf. The differences are returned on output.
   newt uses a globally convergent Newton's method to adjust the values of v until the differences are zero. A user supplied 
   function supplies derivative information to the ODE integrator.
   */
  public double[] funk(final double[] v) {
    double[] v2 = new double[nvar-n2];
    for(int i=0;i<v2.length;i++)
      v2[i] = v[n2+i];
    h1=(x2-x1)/100.0;
    y=buildVector(load1(x1,v));
    Output out = new Output();
    StepperDopr853 s1 = new StepperDopr853();
    Odeint integ1 = new Odeint(y,x1,xf,atol,rtol,h1,hmin,out,d,s1);
    integ1.integrate();
    f1=buildVector(score(xf,y));
    y=buildVector(load2(x2,v2));
    
    StepperDopr853 s2 = new StepperDopr853();
    Odeint integ2 = new Odeint(y,x2,xf,atol,rtol,h1,hmin,out,d,s2);
    integ2.integrate();
    f2=buildVector(score(xf,y));
    for (int i=0;i<nvar;i++) f1[i] -= f2[i];
    return f1;
  }
  
  public abstract double[] load1(final double x, final double[] v);
  
  public abstract double[] load2(final double x, final double[] v);
  
  public abstract double[] score(final double x, final double[] v);
}
