package com.nr.sp;

// Power spectral estimation using overlapping data segments. The user sends non-overlapping 
// segments of length M, which are processed in pairs of 2M, with overlap
public class Spectolap extends Spectreg {
  
  int first;
  double[] fullseg;

  public Spectolap(final int em) {
    super(em);
    first =1;
    fullseg = new double[2*em];
  }

  // Process a data segment of length M using the window function wf
  public void adddataseg(final double[] data, final WindowFun wf) {
    int i;
    if (data.length != m) throw new IllegalArgumentException("wrong size data segment");
    if (first!=0) {   // first segment is just stored
      for (i=0;i<m;i++) fullseg[i+m] = data [i];
      first = 0;
    } else {    // subsequent segments are processed 
      for (i=0;i<m;i++) {
        fullseg[i] = fullseg[i+m];
        fullseg[i+m] = data [i];
      }
      super.adddataseg(fullseg,wf);   // Base class method, the data length is 2M
    }
  }

  // Process a long vector of data as overlapping segments of length 2M
  public void addlongdata(final double[] data, final WindowFun wf) {
    int i, k, noff, nt=data.length, nk=(nt-1)/m;
    double del = nk > 1 ? (nt-m2)/(nk-1.) : 0.;   // Target separation
    if (nt < m2) throw new IllegalArgumentException("data length too short");
    for (k=0;k<nk;k++) {   // Process nk overlapping segments
      noff = (int)(k*del+0.5);   // offset is nearest integer
      for (i=0;i<m2;i++) fullseg[i] = data[noff+i];
      super.adddataseg(fullseg,wf);
    }
  }
}
