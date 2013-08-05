package com.nr.sp;

import static com.nr.NRUtil.*;
import com.nr.fft.FFT;

// Object for accumulating power spectrum estimates from one or more segments of data
public class Spectreg {
  int m,m2,nsum;
  double[] specsum, wksp;

  // Constructor. Sets M, such that data segments will have length 2M, and the spectrum will be estimated at M+1 frequencies
  public Spectreg(final int em){
    m = em;
    m2 =2*m;   // data segment length
    nsum=0;
    specsum = new double[m+1];
    wksp = new double[m2];
    
    if((m & (m-1))!=0) throw new IllegalArgumentException("m must be power of 2");
  }

  // process a data segment of length 2M using the window function wf
  public void adddataseg(final double[] data, final WindowFun wf) {
    int i;
    double w,fac,sumw = 0.;
    if (data.length != m2) throw new IllegalArgumentException("wrong size data segment");
       
    // multiply data with the window function
    for (i=0;i<m2;i++) {  // across signal length
      w = wf.window(i,m2);
      wksp[i] = w*data[i];
      sumw += SQR(w);
    }
    fac = 2./(sumw*m2);
    FFT.realft(wksp,1);
    specsum[0] += 0.5*fac*SQR(wksp[0]);
    for (i=1;i<m;i++) specsum[i] += fac*(SQR(wksp[2*i])+SQR(wksp[2*i+1]));
    specsum[m] += 0.5*fac*SQR(wksp[1]);
    nsum++;
  }

  // return power spectrum estimates as a vector
  public double[] spectrum() {
    double[] spec = new double[m+1];
    if (nsum == 0) throw new IllegalArgumentException("no data yet");
    for (int i=0;i<=m;i++) spec[i] = specsum[i]/nsum;
    return spec;
  }

  // return vector of frequencies (in units of 1/Δ) at which estimates are made
  public  double[] frequencies() {
    double[] freq = new double[m+1];
    for (int i=0;i<=m;i++) freq[i] = i*0.5/m;
    return freq;
  }
}
