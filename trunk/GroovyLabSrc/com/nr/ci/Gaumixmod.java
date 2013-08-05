package com.nr.ci;

import static com.nr.NRUtil.*;
import static java.lang.Math.*;
import com.nr.la.Cholesky;

/*
 Solve a Gaussian mixture model from a set of data points and initial guesses of k means
 */
public class Gaumixmod {
  public int nn, kk, mm;  // numbers of data points, components, and dimensions
  public double[][] data, means, resp;  // local copies of x_n's \mu_k's and the p_{nk}'s
  public double[] frac, lndets;  // P(k)'s and log \Sigma_{k}'s
  public double[][][] sig;
  public double loglike;
  
  // arguments are the data points (as rows in a matrix) and initial 
  // guesses for the means (also as rows in a matrix)
  public Gaumixmod(final double[][] ddata, final double[][] mmeans) {
    int mmstat = ddata[0].length;
    nn = ddata.length;
    kk = mmeans.length;
    mm = mmstat;
    data = buildMatrix(ddata);
    means = buildMatrix(mmeans);
    resp = new double[nn][kk];
    frac = new double[kk];
    lndets = new double[kk];
    sig = new double[kk][mmstat][mmstat];
    
    int i,j,k;
    for (k=0;k<kk;k++) {
      frac[k] = 1./kk;
      for (i=0;i<mm;i++) {
        for (j=0;j<mm;j++) sig[k][i][j] = 0.;
        sig[k][i][i] = 1.0e-10;
      }
    }
    /*
     Perform an initial E-step and M-step  User is responsible for calling additional
     steps until convergence is obtained.
    */ 
    estep();   
    mstep();
  }
  
  public double estep() {
    int k,m,n;
    double tmp,sum,max,oldloglike;
    double[] u = new double[mm],v = new double[mm];
    oldloglike = loglike;
    for (k=0;k<kk;k++) {  // outer loop for computing the p_{nk}'s
 //  Decompose \sigma_{k} in the outer loop     
      Cholesky choltmp = new Cholesky(sig[k]);  
      lndets[k] = choltmp.logdet();
      for (n=0;n<nn;n++) {
        for (m=0;m<mm;m++) u[m] = data[n][m]-means[k][m];
        choltmp.elsolve(u,v);
        for (sum=0.,m=0; m<mm; m++) sum += SQR(v[m]);
        resp[n][k] = -0.5*(sum + lndets[k]) + log(frac[k]);
      }
    }
    loglike = 0;
    for (n=0;n<nn;n++) {  // Seperate normalization for each n
      max = -99.9e99;   // Log-sum-exp trick begins here
      for (k=0;k<kk;k++) if (resp[n][k] > max) max = resp[n][k];
      for (sum=0.,k=0; k<kk; k++) sum += exp(resp[n][k]-max);
      tmp = max + log(sum);
      for (k=0;k<kk;k++) resp[n][k] = exp(resp[n][k] - tmp);
      loglike +=tmp;
    }
    return loglike - oldloglike;  // when abs of this is small, then we have converged 
  }
  
  public void mstep() {
    int j,n,k,m;
    double wgt,sum;
    for (k=0;k<kk;k++) {
      wgt=0.;
      for (n=0;n<nn;n++) wgt += resp[n][k];
      frac[k] = wgt/nn;
      for (m=0;m<mm;m++) {
        for (sum=0.,n=0; n<nn; n++) sum += resp[n][k]*data[n][m];
        means[k][m] = sum/wgt;
        for (j=0;j<mm;j++) {
          for (sum=0.,n=0; n<nn; n++) {
            sum += resp[n][k]*
              (data[n][m]-means[k][m])*(data[n][j]-means[k][j]);
          }
          sig[k][m][j] = sum/wgt;
        }
      }
    }
  }
}
