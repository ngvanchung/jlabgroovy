package com.nr.sort;

import static java.lang.Math.*;
import static com.nr.NRUtil.*;

public class Sorter {
  private Sorter(){}
  
  public static void sort(final double[] arr){
    sort(arr, -1);
  }
  
  // Sort an array arr[0..n-1] into ascending numerical order using the Quicksort
  // algorithm. arr is replaced on output by its sorted rearrangement. Normally, the optional 
  // argument m sould be omitted, but if it is set to a positive value, then only the first m elements
  // of arr are sorted
  public static void sort(final double[] arr, final int m){
    // Here M is the size of subarrays sorted by straight insertion and NSTACK is the required auxiliary storage
      final int M=7, NSTACK=64;
    int i,ir,j,k,jstack=-1,l=0,n=arr.length;
    double a;
    int[] istack = new int[NSTACK];
    if (m>0) n = min(m,n);   // Usew optional argument
    ir=n-1;
    for (;;) {   // Insertion sort when subarray small enough
      if (ir-l < M) {
        for (j=l+1;j<=ir;j++) {
          a=arr[j];
          for (i=j-1;i>=l;i--) {
            if (arr[i] <= a) break;
            arr[i+1]=arr[i];
          }
          arr[i+1]=a;
        }
        if (jstack < 0) break;
        ir=istack[jstack--];   // Pop stack and begin a new round of partitioning
        l=istack[jstack--];   
      } else {
        k=(l+ir) >> 1;       // Choose median of left, center, and right elements as partitioning element a.
        swap(arr, k,l+1);    // Also rearrange so that a[l] <= a[l+1] <= a[ir]
        if (arr[l] > arr[ir]) {
          swap(arr,l,ir);
        }
        if (arr[l+1] > arr[ir]) {
          swap(arr, l+1, ir);
        }
        if (arr[l] > arr[l+1]) {
          swap(arr,l,l+1);
        }
        i=l+1;      // Initialize pointers for partitioning
        j=ir;
        a=arr[l+1];     // Partiiitoning element.
        for (;;) {     // Beginning of innermost loop
          do i++; while (arr[i] < a);    // scan up to find element > a
          do j--; while (arr[j] > a);   // scan down to find element < a
          if (j < i) break;    // pointers crossed. Partitioning complete.
          swap(arr,i,j);   // exchange elements
        }    // end of innermost loop
        arr[l+1]=arr[j];    // insert partitioning element
        arr[j]=a;
        jstack += 2;
          // push pointers to larger subarray on stack; process smaller subarray immediately
        if (jstack >= NSTACK) throw new IllegalArgumentException("NSTACK too small in sort.");
        if (ir-i+1 >= j-l) {
          istack[jstack]=ir;
          istack[jstack-1]=i;
          ir=j-1;
        } else {
          istack[jstack]=j-1;
          istack[jstack-1]=l;
          l=i;
        }
      }
    }
  }
  
  
  // Sort an array arr[0..n-1] into ascending order using Quicksort, while making the corresponding
  // rearrangement of the array brr[0..n-1]
  public static void sort2(final double[] arr, final double[] brr) {
    final int M=7,NSTACK=64;
    int i,ir,j,k,jstack=-1,l=0,n=arr.length;
    double a;
    double b;
    int[] istack = new int[NSTACK];
    ir=n-1;
    for (;;) {   // Insertion sort when subaray small enough
      if (ir-l < M) {
        for (j=l+1;j<=ir;j++) {
          a=arr[j];
          b=brr[j];
          for (i=j-1;i>=l;i--) {
            if (arr[i] <= a) break;
            arr[i+1]=arr[i];
            brr[i+1]=brr[i];
          }
          arr[i+1]=a;
          brr[i+1]=b;
        }
        if (jstack < 0) break;
        ir=istack[jstack--];   // Pop stack and begin a new round of partitioning
        l=istack[jstack--];
      } else {
        k=(l+ir) >> 1;    // Choose median of left, center, and right elements as partitioning element a.
        swap(arr,k,l+1);  // Also rearrange so that a[l] <= a[l+1] <= a[ir]
        swap(brr,k,l+1);
        if (arr[l] > arr[ir]) {
          swap(arr,l,ir);
          swap(brr,l,ir);
        }
        if (arr[l+1] > arr[ir]) {
          swap(arr,l+1,ir);
          swap(brr,l+1,ir);
        }
        if (arr[l] > arr[l+1]) {
          swap(arr,l,l+1);
          swap(brr,l,l+1);
        }
        i=l+1;   // Initialize pointers for partitioning
        j=ir;
        a=arr[l+1];   // Partitioning element
        b=brr[l+1];
        for (;;) {   // Beginning of innermost loop
          do i++; while (arr[i] < a);     // Scan up to find element > a
          do j--; while (arr[j] > a);   // Scan down to find element < a
          if (j < i) break;    // Pointers crossed. Partitioning complete
          swap(arr,i,j);   // Exchange elements of both arrays
          swap(brr,i,j);
        }   // end of innermost loop
        arr[l+1]=arr[j];    // insert partitioning element in both arrays
        arr[j]=a;
        brr[l+1]=brr[j];
        brr[j]=b;
        jstack += 2;
        // push pointers to larger subarray on stack; process smaller subarray subarray immediately
        if (jstack >= NSTACK) throw new IllegalArgumentException("NSTACK too small in sort2.");
        if (ir-i+1 >= j-l) {
          istack[jstack]=ir;
          istack[jstack-1]=i;
          ir=j-1;
        } else {
          istack[jstack]=j-1;
          istack[jstack-1]=l;
          l=i;
        }
      }
    }
  }
  
  
  // sort an array ra[0..n-1] into ascending numerical order using the Heapsort algorithm.
  // ra is replaced on output by its sorted rearrangement
  public static void hpsort(final double[] ra){
    int i,n=ra.length;
    for (i=n/2-1; i>=0; i--)
          // the index i, which here determines the "left" range of the sift-down, i.e., the element
         // to be sifted down, is decremented from n/2-1 dwn to 0 during the "hiring" (heap creation) phase
      sift_down(ra,i,n-1);
    for (i=n-1; i>0; i--) {
        // Here the "right" range of the sift-down is decremented from n-2 down to 0 during the 
        // "retirement-and-promotion" (heap selection) phase
      swap(ra,0,i);   // clear a space a the end of the array, and retire
      sift_down(ra,0,i-1);   // the top of te heap into it
    }
  }
  
  // carry out the sift-down on element ra(l) to maintain the heap structure. l and r
  // determine the "left" and "right" range of the sift-down
  private static void sift_down(final double[]ra, final int l, final int r){
    int j,jold;
    double a;
    a=ra[l];
    jold=l;
    j=2*l+1;
    while (j <= r) {
      if (j < r && ra[j] < ra[j+1]) j++;
      if (a >= ra[j]) break;
      ra[jold]=ra[j];
      jold=j;
      j=2*j+1;
    }
    ra[jold]=a;
  }

  // sort an array arr[0..n-1] into asending numerical order, by straight insertion. arr is replaced 
  // on output by its sorted rearrangement
  // Straight insertion is an O(N^2) routine, and should be used only for small N 
  public static void piksrt(double[] arr){
    int i,j,n=arr.length;
    double a;
    for (j=1;j<n;j++) {
      a=arr[j];
      i=j;
      while (i > 0 && arr[i-1] > a) {
        arr[i]=arr[i-1];
        i--;
      }
      arr[i]=a;
    }
  }

  // sort an array arr[0..n-1] into ascending numerical order, y straight insertion, while making the
  // corresponding rearrangement of the array brr[0..n-1]
  public static void piksr2(double[] arr, double[] brr) {
    int i,j,n=arr.length;
    double a;
    double b;
    for (j=1;j<n;j++) {
      a=arr[j];
      b=brr[j];
      i=j;
      while (i > 0 && arr[i-1] > a) {
        arr[i]=arr[i-1];
        brr[i]=brr[i-1];
        i--;
      }
      arr[i]=a;
      brr[i]=b;
    }
  }
  
  
  // Given k in [0..n-1] returns an array value from arr[0..n-1] such that k array values are
  // less than or equal to the one returned. The input array will be rearranged to have this 
  // value in location arr[k], with all smaller elements moved to arr[0..k-1] (in arbitrary order)
  // and all larger elements in arr[k+1..n-1] (also in arbitrary order)
  public static double select(final int k, final double[] arr) {
    int i,ir,j,l,mid,n=arr.length;
    double a;
    l=0;
    ir=n-1;
    for (;;) {  
      if (ir <= l+1) {    //  Active partition contains 1 or 2 elements 
        if (ir == l+1 && arr[ir] < arr[l]){    // Case of 2 elements
          swap(arr,l,ir);
        }
        return arr[k];
      } else {
        mid=(l+ir) >> 1;       // Choose median of left, center, and right elements as partitioning element a.
        swap(arr,mid,l+1);   // Also, rearrange so that arr[l] <= arr[l+1], arr[ir] >= arr[l+1]
        if (arr[l] > arr[ir]){
          swap(arr, l, ir);
        }
        if (arr[l+1] > arr[ir]){
          swap(arr, l+1, ir);
        }
        if (arr[l] > arr[l+1]){
          swap(arr,l,l+1);
        }
        i=l+1;     // Initialize pointers for partitioning.
        j=ir;
        a=arr[l+1];     // Partitioning element.
        for (;;) {   // Beginning of innermost loop.
          do i++; while (arr[i] < a);   // Scan up to find element > a
          do j--; while (arr[j] > a);   // Scan down to find element < a
          if (j < i) break;    // Pointers crossed. Partitioning complete
          swap(arr,i,j);
        }   // end of innermost loop
        arr[l+1]=arr[j];    // insert partitioning element
        arr[j]=a;
        if (j >= k) ir=j-1;   // keep active the partition that contains the kth element
        if (j <= k) l=i;
      }
    }
  }
  
  public static void shell(final double[]a){
    shell(a, -1);
  }
  
  // Sort an array a[0..n-1] into ascending numerical order by Shell's method (diminishing increment sort)
  // a is replaced on output by its sorted rearrangement. Normally, the optional argument m should be ommited,
  // but if it is set to a positive value, then only the first m elements of a are sorted
  public static void shell(final double[]a, final int m){
    int i,j,inc,n=a.length;
    double v;
    if (m>0) n = Math.min(m,n);
    inc=1;
    do {
      inc *= 3;
      inc++;
    } while (inc <= n);
    do {
      inc /= 3;
      for (i=inc;i<n;i++) {
        v=a[i];
        j=i;
        while (a[j-inc] > v) {
          a[j]=a[j-inc];
          j -= inc;
          if (j < inc) break;
        }
        a[j]=v;
      }
    } while (inc > 1);
  }
  
  
  /**
   * select Mth largest in place
   * 
   */
  public static double selip(final int k, final double[] arr) {
    final int M=64;
    final double BIG=.99e99;
    int i,j,jl,jm,ju,kk,mm,nlo,nxtmm,n=arr.length;
    double ahi,alo,sum;
    int[] isel = new int[M+2];
    double[] sel = new double[M+2];
    if (k < 0 || k > n-1) throw new IllegalArgumentException("bad input to selip");
    kk=k;
    ahi=BIG;
    alo = -BIG;
    for (;;) {
      mm=nlo=0;
      sum=0.0;
      nxtmm=M+1;
      for (i=0;i<n;i++) {
        if (arr[i] >= alo && arr[i] <= ahi) {
          mm++;
          if (arr[i] == alo) nlo++;
          if (mm <= M) sel[mm-1]=arr[i];
          else if (mm == nxtmm) {
            nxtmm=mm+mm/M;
            sel[(i+2+mm+kk) % M]=arr[i];
          }
          sum += arr[i];
        }
      }
      if (kk < nlo) {
        return alo;
      }
      else if (mm < M+1) {
        Sorter.shell(sel,mm);
        ahi = sel[kk];
        return ahi;
      }
      sel[M]=sum/mm;
      Sorter.shell(sel,M+1);
      sel[M+1]=ahi;
      for (j=0;j<M+2;j++) isel[j]=0;
      for (i=0;i<n;i++) {
        if (arr[i] >= alo && arr[i] <= ahi) {
          jl=0;
          ju=M+2;
          while (ju-jl > 1) {
            jm=(ju+jl)/2;
            if (arr[i] >= sel[jm-1]) jl=jm;
            else ju=jm;
          }
          isel[ju-1]++;
        }
      }
      j=0;
      while (kk >= isel[j]) {
        alo=sel[j];
        kk -= isel[j++];
      }
      ahi=sel[j];
    }
  }
}
