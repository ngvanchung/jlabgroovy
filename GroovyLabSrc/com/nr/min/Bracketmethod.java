package com.nr.min;

import static com.nr.NRUtil.*;
import static java.lang.Math.*;
import com.nr.UniVarRealValueFun;

/**
 * Base class for one-dimensional minimization routines. Provides a routine to
 * bracket a minimum and several utility functions.
 * 
 */
public class Bracketmethod {
  public double ax, bx, cx, fa, fb, fc;
  
  public Bracketmethod(){
    
  }
  
  /**
   * Given a function or functor func, and given distinct initial points ax and
   * bx, this routine searches in the downhill direction (defined by the
   * function as evaluated at the initial points) and returns new points ax, bx,
   * cx that bracket a minimum of the function. Also returned are the function
   * values at the three points, fa, fb, and fc.
   * 
   * @param a
   * @param b
   * @param func
   */
  public void bracket(final double a, final double b, final UniVarRealValueFun func) {
    // here GOLD is the default ratio by which successive intervals are magnified and GLIMIT
    // is the maximum magnification allowed for a parabolic-fit step  
    final double GOLD=1.618034, GLIMIT=100.0, TINY=1.0e-20;
    ax=a; bx=b;
    double fu;
    fa=func.funk(ax);
    fb=func.funk(bx);
    if (fb > fa) {  // Switch roles of a and b so that we can go downhill in the direction from a to b
      //SWAP(ax,bx);
      double swap=ax; ax=bx; bx=swap;
      //SWAP(fb,fa);
      swap=fb; fb=fa; fa=swap;
    }
    cx=bx+GOLD*(bx-ax);   // first guess for c
    fc=func.funk(cx);
    while (fb > fc) {    // keep returning here until we bracket
      double r=(bx-ax)*(fb-fc);     // compute u by parabolic extrapolation from a, b, c.
      double q=(bx-cx)*(fb-fa);   // TINY is used to prevent any possible division by zero
      double u=bx-((bx-cx)*q-(bx-ax)*r)/
        (2.0*SIGN(max(abs(q-r),TINY),q-r));
      double ulim=bx+GLIMIT*(cx-bx);
      // we won't go farther than this. Test various possibilities:
      if ((bx-u)*(u-cx) > 0.0) {   // Parabolic u is between b and c: try it
        fu=func.funk(u);
        if (fu < fc) {    // Got a minimum between b and c.
          ax=bx;
          bx=u;
          fa=fb;
          fb=fu;
          return;
        } else if (fu > fb) {   // Got a minimum between a and u
          cx=u;
          fc=fu;
          return;
        }
        u=cx+GOLD*(cx-bx);   // Parabolic fit was no use. Use default magnification.
        fu=func.funk(u);
      } else if ((cx-u)*(u-ulim) > 0.0) {   // Parabolic fit is between c and its allowed limit.
        fu=func.funk(u);
        if (fu < fc) {
          //shft3(bx,cx,u,u+GOLD*(u-cx));
          double dum = u+GOLD*(u-cx);
          bx=cx; cx=u; u=dum;
          //shft3(fb,fc,fu,func.funk(u));
          fb = fc; fc = fu; fu=func.funk(u);
        }
      } else if ((u-ulim)*(ulim-cx) >= 0.0) {   // Limit parabolic u to maximum allowed value
        u=ulim;
        fu=func.funk(u);
      } else {    // Reject parabolic u, use default magnification
        u=cx+GOLD*(cx-bx);
        fu=func.funk(u);
      }
      //shft3(ax,bx,cx,u);
      ax=bx; bx=cx; cx=u;   // eliminate oldest point and continue
      //shft3(fa,fb,fc,fu);
      fa=fb; fb=fc; fc=fu;
    }
  }
  
  /*
  public static void shft2(doubleW a, doubleW b, final double c) {
    a.val=b.val;
    b.val=c;
  }
  
  public static void shft3(doubleW a, doubleW b, doubleW c, final double d)  {
    a.val=b.val;
    b.val=c.val;
    c.val=d;
  }
  
  public static void mov3(doubleW a, doubleW b, doubleW c, final double d, final double e, final double f) {
    a.val=d; b.val=e; c.val=f;
  }
  */
}
