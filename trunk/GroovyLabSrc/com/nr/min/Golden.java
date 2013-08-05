package com.nr.min;

import static java.lang.Math.*;
import com.nr.UniVarRealValueFun;

/**
 * Golden Section Search in One Dimension
 * 
 * Copyright (C) Numerical Recipes Software 1986-2007
 * Java translation Copyright (C) Huang Wen Hui 2012
 *
 * @author hwh
 *
 */

// Golden search for minimum
/*
 Given a function f, and given a bracketing triplet of abscissas ax, bx, cx
 (such that bx is between ax and cx, and f(bx) is less than both f(ax)
 and f(cx)), this routine performs a golden section search for the minimum,
 isolating it to a fractional precision of about tol. The abscissa of the minimum 
 is returned as xmin, and the function value at the minimum as min, the returned
 function value.
 */
public class Golden extends Bracketmethod{
  double xmin,fmin;
  final double tol;
  
  public Golden(){
    this(3.0e-8);
  }
  
  public Golden(final double toll) {
    tol = toll;
  }

  public double minimize(final UniVarRealValueFun func)  {
    final double R=0.61803399,C=1.0-R;   // the golden ratios
    double x1,x2;
    double x0=ax;   // at any given time we will keep track of four points, x0, x1, x2, x3
    double x3=cx;
    if (abs(cx-bx) > abs(bx-ax)) {   // Make x0 to x1 the smaller segment
      x1=bx;
      x2=bx+C*(cx-bx);    // and fill in the new point to be tried
    } else {
      x2=bx;
      x1=bx-C*(bx-ax);
    }
    double f1=func.funk(x1);   // the initial function evaluations. Note that 
    double f2=func.funk(x2);    // we never need to evaluate the function
    while (abs(x3-x0) > tol*(abs(x1)+abs(x2))) {    //// at the original endpoints........
      if (f2 < f1) {     //  one possible outcome,
        //shft3(x0,x1,x2,R*x2+C*x3);  // its housekeeping,
        double dum = R*x2+C*x3;
        x0=x1; x1=x2; x2=dum;
        // shft2(f1,f2,func.funk(x2));   // and a new function evaluation
        f1=f2; f2 = func.funk(x2);
      } else {     // the other outcome
        //shft3(x3,x2,x1,R*x1+C*x0);
        double dum = R*x1+C*x0;
        x3=x2; x2=x1; x1=dum;
        // shft2(f2,f1,func.funk(x1));    // and its new function evaluation
        f2=f1;f1=func.funk(x1);
      }
    }    // Back to see if we are done.
    if (f1 < f2) {     // We are done. Output the best of the two 
      xmin=x1;    // current valuesss.
      fmin=f1;
    } else {
      xmin=x2;
      fmin=f2;
    }
    return xmin;
  }
}
