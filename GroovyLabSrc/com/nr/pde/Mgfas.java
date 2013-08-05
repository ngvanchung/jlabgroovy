package com.nr.pde;

import static com.nr.NRUtil.*;
import static java.lang.Math.*;

import org.netlib.util.doubleW;

/*"
 Full multigrid algorithm for FAS solution of nonlinear elliptic equation, here equation (20.6.44)
 on a square domain of side 1, so that h = 1/(n-1)
 */
public class Mgfas {
  int n,ng;
  double[][] uj,uj1;
  //NRvector<NRmatrix<double> *> rho;
  double[][][] rho;   // Vector of pointers to ρ on each level
  
  public double[][] u;

  /**
   * 
   * @param u
   * @param maxcyc
   */
  /*
   On input u[0..n-1][0..n-1] contains the right-hand side ρ, while on output it returns the
   solution. The dimension n must be of the form 2^j+1 for some integer j. (j is actually the number 
   of grid levels used in the solution, called ng below). maxcyc is the maximum number of V-cycles to be used at each level
   */
  public Mgfas(final double[][] u, final int maxcyc) {
    this.u = u;
    n = u.length;
    ng = 0;
    int nn=n;
    while ((nn >>>= 1)!=0) ng++;
    if ((n-1) != (1 << ng))
      throw new IllegalArgumentException("n-1 must be a power of 2 in mgfas.");
    nn=n;
    int ngrid=ng-1;
    rho = new double[ng][][];   // Allocate storage for r, h, s on grid ng-1
    //rho.resize(ng);                  // and fill it with the input r, h, s
    rho[ngrid]=new double[nn][nn];   // Similarly allocate storage and fill r, h, s by
    copyAssign(rho[ngrid],u);           // restriction on all coarse grids
    while (nn > 3) {
      nn=nn/2+1;
      rho[--ngrid]=new double[nn][nn];
      rstrct(rho[ngrid],rho[ngrid+1]);
    }
    nn=3;
    uj=new double[nn][nn];
    slvsm2(uj,rho[0]);    // Initial solution on coarsest grid
    for(int j=1;j<ng;j++) {  // nested iteration loop
      nn=2*nn-1;
      uj1=uj;
      uj=new double[nn][nn];
      double[][] temp = new double[nn][nn];
      interp(uj,uj1);  // interpolate from grid j-1 to next finer grid j
      for (int jcycle=0;jcycle<maxcyc;jcycle++) {        // V-cycle loop
        doubleW trerr= new doubleW(1.0);   // R.h.s is dummy
        mg(j,uj,temp,rho,trerr);
        lop(temp,uj);   // Form residual ||d_h||
        matsub(temp,rho[j],temp);
        double res=anorm2(temp);
        if (res < trerr.val) break;    // No more V-cycles needed if residual small enough
      }
    }
    this.u = uj;   // return solution in u
  }
  
  
  /*
   Matrix addition: Adds a[0..n-1][0..n-1] to b[0..n-1][0..n-1] and returns reslt in c[0..n-1][0..n-1]
   */
  public static void matadd(final double[][] a, final double[][] b, final double[][] c) {
    int n=a.length;
    for (int j=0;j<n;j++)
      for (int i=0;i<n;i++)
        c[i][j]=a[i][j]+b[i][j];
  }
  
  /*
   Matrix subtraction: Subtracts b[0..n-1][0..n-1] from a[0..n-1][0..n-1] and returns result in c[0..n-1][0..n-1]
   */
  public static void matsub(final double[][] a, final double[][] b, final double[][] c) {
    int n=a.length;
    for (int j=0;j<n;j++)
      for (int i=0;i<n;i++)
        c[i][j]=a[i][j]-b[i][j];
  }
  
/* 
   Solution of equation (20.6.44) on the coarsest grid, where h = 1/2. The right-hand side is 
   input in rhs[0..2][0..2] and the solution is returned in u[0..2][0..2]
   */
  public static void slvsm2(final double[][] u, final double[][] rhs){
    double h=0.5;
    for (int i=0;i<3;i++)
      for (int j=0;j<3;j++)
        u[i][j]=0.0;
    double fact=2.0/(h*h);
    double disc=sqrt(fact*fact+rhs[1][1]);
    u[1][1]= -rhs[1][1]/(fact+disc);
  }
  
  /*
   Red-black Gauss-Seidel relaxation for equation (20.6.44). The current value of the solution 
   u[0..n-1][0..n-1] is updated, using the right-hand side function rhs[0..n-1][0..n-1]
   */
  public static void relax2(final double[][] u, final double[][] rhs) {
    int n=u.length;
    int jsw=1;
    double h=1.0/(n-1);
    double h2i=1.0/(h*h);
    double foh2 = -4.0*h2i;
    for (int ipass=0;ipass<2;ipass++,jsw=3-jsw) {
      int isw=jsw;
      for (int j=1;j<n-1;j++,isw=3-isw) {
        for (int i=isw;i<n-1;i+=2) {
          double res=h2i*(u[i+1][j]+u[i-1][j]+u[i][j+1]+u[i][j-1]-
            4.0*u[i][j])+u[i][j]*u[i][j]-rhs[i][j];
          u[i][j] -= res/(foh2+2.0*u[i][j]);
        }
      }
    }
  }
  
  /*
   Half-weighting restriction. If nc is the coarse-grid dimension, the fine-grid solution is input 
   in uf[0..2*nc-2][0..2*nc-2]. The coarse-grid solution obtained by restriction is returned in uc[0..nc-1][0..nc-1]
   */
  public static void rstrct(final double[][] uc, final double[][] uf) {
    int nc=uc.length;
    int ncc=2*nc-2;
    for (int jf=2,jc=1;jc<nc-1;jc++,jf+=2) {   // Interior points
      for (int iif=2,ic=1;ic<nc-1;ic++,iif+=2) {
        uc[ic][jc]=0.5*uf[iif][jf]+0.125*(uf[iif+1][jf]+uf[iif-1][jf]
          +uf[iif][jf+1]+uf[iif][jf-1]);
      }
    }
    for (int jc=0,ic=0;ic<nc;ic++,jc+=2) {   // Boundary points
      uc[ic][0]=uf[jc][0];
      uc[ic][nc-1]=uf[jc][ncc];
    }
    for (int jc=0,ic=0;ic<nc;ic++,jc+=2) {
      uc[0][ic]=uf[0][jc];
      uc[nc-1][ic]=uf[ncc][jc];
    }
  }
  
  /*
   Given u[0..n-1][0..n-1], returns \phi_h(u_h) for eqn. (20.6.44) in ot[0..n-1][0..n-1]
   */
  public static void lop(final double[][] out, final double[][] u) {
    int n=u.length;
    double h=1.0/(n-1);
    double h2i=1.0/(h*h);
    for (int j=1;j<n-1;j++)    // Interior points
      for (int i=1;i<n-1;i++)
        out[i][j]=h2i*(u[i+1][j]+u[i-1][j]+u[i][j+1]+u[i][j-1]-
          4.0*u[i][j])+u[i][j]*u[i][j];
    for (int i=0;i<n;i++)    // Boundary points
      out[i][0]=out[i][n-1]=out[0][i]=out[n-1][i]=0.0;
  }
  
  /*
   Coarse-to-fine prolongation by bilinear interpolation. If nf is the fine-grid dimension, the coarse-grid solution is
   input as uc[0..nc-1][0..nc-1], where nc = nf/2+1. The fine-grid solution is returned in f[0..nf-1][0..nf-1]
   */
  public static void interp(final double[][] uf, final double[][] uc){
    int nf=uf.length;
    int nc=nf/2+1;
    for (int jc=0;jc<nc;jc++)   // Do elements that are copies
      for (int ic=0;ic<nc;ic++) uf[2*ic][2*jc]=uc[ic][jc];
    for (int jf=0;jf<nf;jf+=2)   // Do even-numbered columns, interpolating verically
      for (int iif=1;iif<nf-1;iif+=2)
        uf[iif][jf]=0.5*(uf[iif+1][jf]+uf[iif-1][jf]);
    for (int jf=1;jf<nf-1;jf+=2)   // Do odd-numbered columns, interpolating horizontally
      for (int iif=0;iif<nf;iif++)
        uf[iif][jf]=0.5*(uf[iif][jf+1]+uf[iif][jf-1]);
  }
  
  // Returns the Euclidean norm of the matrix a[0..n-1][0..n-1]
  public static double anorm2(final double[][] a) {
    double sum=0.0;
    int n=a.length;
    for (int j=0;j<n;j++)
      for (int i=0;i<n;i++)
        sum += a[i][j]*a[i][j];
    return sqrt(sum)/n;
  }
  
  /*
   Recursive multigrid iteration. On input, j is the current level and u is the current value of the solution.
   For the first call on a given level, the right-hand side is zero, and the argument rhs is dummy. 
   This is signaled by inputting trerr positive. Subsequent recursive calls supply
   a nonzero rhs as in equation (20.6.33). Tis is signaled by inputting trerr negative. rho is the vector 
   of pointers to ρ on each level. On output u contains the improved solution at the current level. When the
   first call on a given level is made, the relative truncation error is returned in trerr
   */
  public static void mg(final int j, final double[][] u, final double[][] rhs,
    final double[][][]rho, final doubleW trerr) {
    final int NPRE=1,NPOST=1;
    // Number of relaxation sweeps before and after the coarse-grid correction is computed
    final double ALPHA=0.33;   // Relates the estimated truncation error to the norm of the residual
    doubleW dum= new doubleW(-1.0);
    int nf=u.length;
    int nc=(nf+1)/2;
    double[][] temp = new double[nf][nf];
    if (j == 0) {   // Bottom of V:  Solve on coarsest grid.
      matadd(rhs,rho[j],temp);
      slvsm2(u,temp);
    } else {   // On downward stoke of the V
      double[][] v= new double[nc][nc],ut=new double[nc][nc],tau=new double[nc][nc],tempc=new double[nc][nc];
      for (int jpre=0;jpre<NPRE;jpre++) {   // Pre-smoothing
        if (trerr.val < 0.0) {  
          matadd(rhs,rho[j],temp);
          relax2(u,temp);
        }
        else
          relax2(u,rho[j]);
      }
      rstrct(ut,u);
      copyAssign(v, ut);
      lop(tau,ut);
      lop(temp,u);
      if (trerr.val < 0.0)
        matsub(temp,rhs,temp);
      rstrct(tempc,temp);
      matsub(tau,tempc,tau);
      if (trerr.val > 0.0)
        trerr.val=ALPHA*anorm2(tau);
      mg(j-1,v,tau,rho,dum);
      matsub(v,ut,tempc);
      interp(temp,tempc);
      matadd(u,temp,u);
      for (int jpost=0;jpost<NPOST;jpost++) {
        if (trerr.val < 0.0) {
          matadd(rhs,rho[j],temp);
          relax2(u,temp);
        }
        else
          relax2(u,rho[j]);
      }
    }
  }
}
