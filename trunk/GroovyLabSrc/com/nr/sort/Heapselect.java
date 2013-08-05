package com.nr.sort;

import static com.nr.NRUtil.*;
import static java.lang.Math.min;

// Class for tracking the m largest values seen thus far in a stream of values
public class Heapselect {
  int m,n,srtd;
  double[] heap;

  // the argument is the number of largest values to track
  public Heapselect(int mm) {
    m =mm;
    n=0;
    srtd=0;
    heap = buildVector(mm, 1.e99);
  }

  // assiiiiiilate a new value from the stream
  public void add(double val) {
    int j,k;
    if (n<m) {   // Heap not yet filled
      heap[n++] = val;  
      if (n==m) Sorter.sort(heap);   // Create initial heap by overkill!
    } else {
      if (val > heap[0]) {   // Put it on the heap?
        heap[0]=val;
        for (j=0;;) {   // Sift down
          k=(j << 1) + 1;
          if (k > m-1) break;
          if (k != (m-1) && heap[k] > heap[k+1]) k++;
          if (heap[j] <= heap[k]) break;
          swap(heap,k,j);
          j=k;
        }
      }
      n++;
    }
    srtd = 0;   // Mark heap as "unsorted"
  }

  // Return the kth largest value seen so far. k = 0 returns the largest value seen, k = 1 the second largest, ...,
  // k = m-1 the last position tracked. Also, k must be less than the number of previous values assimilated
  public double report(int k) {
    int mm = min(n,m);
    if (k > mm-1) throw new IllegalArgumentException("Heapselect k too big");
    if (k == m-1) return heap[0];    // always free, since top of heap
    if (srtd==0) { Sorter.sort(heap); srtd = 1; }   // otherwise, need to sort the heap
    return heap[mm-1-k];
  }
}
