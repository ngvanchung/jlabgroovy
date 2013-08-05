package com.nr.sp;

import org.netlib.util.doubleW;

import com.nr.fft.*;
import static java.lang.Math.*;
import com.nr.UniVarRealValueFun;
import com.nr.interp.Poly_interp;

/**
 * Example program illustrating how to use the routine dftcor. 
 * The user supplies an external function func that returns the
 * quantity h(t). The routine then returns S(a,b, cos(wt)*h(t)*dt as
 * cosint and S(a,b, sin(wt)*h(t)*dt as sinint.
 * 
 */
public class DftInt {
  final int M=64,NDFT=1024,MPOL=6;
  /*
   The values of M, NDFT, and MPOL are merely illustrative and should be optimized for your
   particular application. M is the number of subintervals, NDFT is the length of the FFT ( a power 
   of 2), and MPOL is the degree of polynomial interpolation used to obtain the desired frequency from the FFT
   */
  double aold = -1.e30,bold = -1.e30,delta;
  int init=0;
  double[] data = new double[NDFT];
  double[] endpts = new double[8];
  UniVarRealValueFun funcold;
  

  public void dftint(final UniVarRealValueFun func, final double a, final double b, final double w,
      final doubleW cosint, final doubleW sinint) {
    final double TWOPI=6.283185307179586476;
    int j,nn;
    double c,cdft,en,s,sdft;
    doubleW corfac = new doubleW(0);
    doubleW corim = new doubleW(0);
    doubleW corre = new doubleW(0);
    double[] cpol = new double[MPOL];
    double[] spol = new double[MPOL];
    double[] xpol = new double[MPOL];
    if (init != 1 || a != aold || b != bold || func != funcold) {
        // Do we need to initialize, or is only ω changed?
      init=1;
      aold=a;
      bold=b;
      funcold=func;
      delta=(b-a)/M;
      for (j=0;j<M+1;j++)   //    Load the function values into the data array.
        data[j]=func.funk(a+j*delta);
      for (j=M+1;j<NDFT;j++)   // Zero-pad the rest of the data array
        data[j]=0.0;
      for (j=0;j<4;j++) {    // Load the endpoints.
        endpts[j]=data[j];
        endpts[j+4]=data[M-3+j];
      }
      FFT.realft(data,1);   
      // realft returns the unused value corresponding to ω_{N/2} in data[1]. We actually want
      // this element to contain the imaginary part corresponding to ω_0, which is zero
      data[1]=0.0;
    }
    // Now interpolate on the DFT result for the desired frequency. If the frequency is an ω_n,
    // i.e. the quantity en is an integer, then cdft = data[2*en-2], sdft = ata[2*en-1], 
    // and you could omit the interpolation
    en=w*delta*NDFT/TWOPI+1.0;
    nn=min(max((int)(en-0.5*MPOL+1.0),1),NDFT/2-MPOL+1);  // Leftmost point for the interpolation
    for (j=0;j<MPOL;j++,nn++) {
      cpol[j]=data[2*nn-2];
      spol[j]=data[2*nn-1];
      xpol[j]=nn;
    }
    cdft = new Poly_interp(xpol,cpol,MPOL).interp(en);
    sdft = new Poly_interp(xpol,spol,MPOL).interp(en);
    Fourier.dftcor(w,delta,a,b,endpts,corre,corim,corfac);    // Now get the endpoint correction and the multiplicative
    cdft *= corfac.val;                                                         // factor W(θ)
    sdft *= corfac.val;
    cdft += corre.val;
    sdft += corim.val;
    c=delta*cos(w*a);       // Finally multiply by Δ and exp(ιωα)
    s=delta*sin(w*a);
    cosint.val=c*cdft-s*sdft;
    sinint.val=s*cdft+c*sdft;
  }
}
