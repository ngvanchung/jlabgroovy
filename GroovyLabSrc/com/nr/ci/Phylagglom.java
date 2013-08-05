package com.nr.ci;

import static com.nr.NRUtil.*;

import java.io.*;

import org.netlib.util.doubleW;

// Abstract base class for constructing an agglomerative phylogenetic ree
public abstract class Phylagglom {
  public int n, root, fsroot;  // No. of data points, root node, forced root
  public double seqmax, depmax;  // max. values of seq, dep over the tree
  public Phylagglomnode[] t;   // the tree
  public abstract void premin(double[][] d, int[] nextp); //  function called before minimum search  
  public abstract double dminfn(double[][] d, int i, int j); //  function to be minimized  
  public abstract double dbranchfn(double[][] d, int i, int j); //  branch length, node i to mother (j is sister)  
  public abstract double dnewfn(double[][] d, int k, int i, int j, int ni, int nj);  // distance function for newly constructed nodes
  public abstract void drootbranchfn(double[][] d, int i, int j, int ni, int nj, doubleW bi, doubleW bj);    // sets branch lengths to the final root node

  public Phylagglom(final double[][] dist){
    this(dist,-1);
  }
  
  public Phylagglom(final double[][] dist, final int fsr) {
    n = dist.length;
    fsroot = fsr;
    t = new Phylagglomnode[2*n-1];
    for(int i=0;i<t.length;i++)
      t[i] = new Phylagglomnode();
  }

  // routine that actually constructs the tree, called by the constructor of a derived class
  public void makethetree(final double[][] dist) {
    int i, j, k, imin=0, jmin=0, node, ntask; 
    double dd, dmin;
    double[][] d = buildMatrix(dist);  // Matrix d is initialized with dist
    int[] tp = new int[n], nextp = new int[n], prevp = new int[n], tasklist = new int[2*n+1];
    double[] tmp = new double[n];
    for (i=0;i<n;i++) {  // initializations on leaf elements
      // nextp and prevp are for looping on the distance matrix even as it becomes sparse
      nextp[i] = i+1;   
      prevp[i] = i-1;   
      
      tp[i] = i;  // tp points from a distance mattttttttix row to a tree element
      t[i].ldau = t[i].rdau = -1;
      t[i].nel = 1;
    }   
    prevp[0] = nextp[n-1] = -1;   // Signifying end of loop
    // ncurr = n;
    for (node = n; node < 2*n-2; node++) {   // Main loop!
      premin(d,nextp);  // Any calculations needed before min finding
      dmin = 9.99e99;
      for (i=0; i>=0; i=nextp[i]) {   // Find i,j pair with min distance
        if (tp[i] == fsroot) continue;
        for (j=nextp[i]; j>=0; j=nextp[j]) {
          if (tp[j] == fsroot) continue;
          if ((dd = dminfn(d,i,j)) < dmin) {
            dmin = dd;
            imin = i; jmin = j;
          }
        }
      }
      i = imin; j = jmin;
        // Now set properties of the parent and children
      t[tp[i]].mo = t[tp[j]].mo = node;
      t[tp[i]].modist = dbranchfn(d,i,j);
      t[tp[j]].modist = dbranchfn(d,j,i);
      t[node].ldau = tp[i];
      t[node].rdau = tp[j];
      t[node].nel = t[tp[i]].nel + t[tp[j]].nel;
      for (k=0; k>=0; k=nextp[k]) {  // Get new-node distances
        tmp[k] = dnewfn(d,k,i,j,t[tp[i]].nel,t[tp[j]].nel);
      }
      for (k=0; k>=0; k=nextp[k]) d[i][k] = d[k][i] = tmp[k];
      // new node replaces child i in dist. matrix, while child gets patched around
      tp[i] = node;  
      if (prevp[j] >= 0) nextp[prevp[j]] = nextp[j];
      if (nextp[j] >= 0) prevp[nextp[j]] = prevp[j];
      
    }  // End of main loop
    
    i = 0; j = nextp[0];  // set properties of the root node
    root = node;
    t[tp[i]].mo = t[tp[j]].mo = t[root].mo = root;
    
    doubleW bi = new doubleW(t[tp[i]].modist);
    doubleW bj = new doubleW(t[tp[j]].modist);
    drootbranchfn(d,i,j,t[tp[i]].nel,t[tp[j]].nel, bi,bj);
      //t[tp[i]].modist,t[tp[j]].modist);
    t[tp[i]].modist = bi.val;
    t[tp[j]].modist = bj.val;
    
    t[root].ldau = tp[i];
    t[root].rdau = tp[j];
    t[root].modist = t[root].dep = 0.;
    t[root].nel = t[tp[i]].nel + t[tp[j]].nel;
    
    // we now traverse the tree computing seq and dep, hints for where to plot nodes in a two-dimensional representation. See Numerical Recipes book text
    ntask = 0;
    seqmax = depmax = 0.;
    tasklist[ntask++] = root;
    while (ntask > 0) {
      i = tasklist[--ntask];
      if (i >= 0) {
        t[i].dep = t[t[i].mo].dep + t[i].modist;
        if (t[i].dep > depmax) depmax = t[i].dep;
        if (t[i].ldau < 0) {
          t[i].seq = seqmax++;
        } else {
          tasklist[ntask++] = -i-1;
          tasklist[ntask++] = t[i].ldau;
          tasklist[ntask++] = t[i].rdau;
        }
      } else {
        i = -i-1;
        t[i].seq = 0.5*(t[t[i].ldau].seq + t[t[i].rdau].seq);
      }
    }
  }
  
  public int comancestor(final int leafa, final int leafb) {
    int i, j;
    for (i = leafa; i != root; i = t[i].mo) {
      for (j = leafb; j != root; j = t[j].mo) if (i == j) break;
      if (i == j) break;
    }
    return i;
  }
  
  /*
   Output a pylogenetic tree p in the "Newick" or "New Hampshire" standard format. Text labels
   for the leaf nodes are input as the rows of str, each a null terminated string. The output file name is specified by filename
   */
  public static void newick(Phylagglom p, char[][] str, String filename) throws IOException {
    java.io.PrintWriter OUT =new java.io.PrintWriter(new FileWriter(filename));
    int i, s, ntask = 0, n = p.n, root = p.root;
    int[] tasklist = new int[2*n+1];
    tasklist[ntask++] = (1 << 16) + root;
    while (ntask-- > 0) {
      s = tasklist[ntask] >> 16;
      i = tasklist[ntask] & 0xffff;
      if (s == 1 || s == 2) {
        tasklist[ntask++] = ((s+2) << 16) + p.t[i].mo;
        if (p.t[i].ldau >= 0) {
          OUT.printf("(");
          tasklist[ntask++] = (2 << 16) + p.t[i].rdau;      
          tasklist[ntask++] = (1 << 16) + p.t[i].ldau;      
        }
        else OUT.printf("%s:%f",new String(str[i], 0, str[i].length-1),p.t[i].modist); 
      }
      else if (s == 3) {if (ntask > 0) OUT.printf(",\n");}
      else if (s == 4) {
        if (i == root) OUT.printf(");\n");
        else OUT.printf("):%f",p.t[i].modist);
      }
    }
    OUT.close();
  }
  
  public static void phyl2ps(String filename, Phylagglom ph, char[][] str, int extend,
    double xl, double xr, double yt, double yb)  throws IOException{
    int i,j;
    double id,jd,xi,yi,xj,yj,seqmax,depmax;
    java.io.PrintWriter OUT =new java.io.PrintWriter(new FileWriter(filename));
    OUT.printf("%%!PS\n/Courier findfont 8 scalefont setfont\n");
    seqmax = ph.seqmax;
    depmax = ph.depmax;
    for (i=0; i<2*(ph.n)-1; i++) {
      j = ph.t[i].mo;
      id = ph.t[i].dep;
      jd = ph.t[j].dep;
      xi = xl + (xr-xl)*id/depmax;
      yi = yt - (yt-yb)*(ph.t[i].seq+0.5)/seqmax;
      xj = xl + (xr-xl)*jd/depmax;
      yj = yt - (yt-yb)*(ph.t[j].seq+0.5)/seqmax;
      OUT.printf("%f %f moveto %f %f lineto %f %f lineto 0 setgray stroke\n",
        xj,yj,xj,yi,xi,yi);
      if (extend!=0) {
        if (i < ph.n) {
          OUT.printf("%f %f moveto %f %f lineto 0.7 setgray stroke\n",
            xi,yi,xr,yi);
          OUT.printf("%f %f moveto (%s (%02d)) 0 setgray show\n",
            xr+3.,yi-2.,new String(str[i], 0, str[i].length-1),i);
        }
      } else {
        if (i < ph.n) OUT.printf("%f %f moveto (%s (%02d)) 0 setgray show\n",
          xi+3.,yi-2.,new String(str[i], 0, str[i].length-1),i);
      }
    }
    OUT.printf("showpage\n\004");
    OUT.close();
  }
}
