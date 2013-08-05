package com.nr.ci;

// Defines what a kernel structure needs to provide
public abstract class Svmgenkernel {
  int m, kcalls;  // No. of data points, counter for kernel calls

  double[][] ker;  // locally stored kernel matrix

  double[] y;  // must provide reference to the y_i's

  double[][] data;  // must provide reference to the x_i's

  public Svmgenkernel(final double[] yy, final double[][] ddata) {
    m = yy.length;
    kcalls = 0;
    ker = new double[m][m];
    y = yy;
    data = ddata;
  }

  // every kernel structure must provide a kernel function that returns the kernel  for arbitrary feature vectors
  public abstract double kernel(final double[] xi, final double[] xj);

  public double kernel(final int i, final double[] xj) {
    return kernel(data[i], xj);
  }

  // every kernel structure's constructor must call fill to fill the ker  matrix 
  public void fill() {
    int i, j;
    for (i = 0; i < m; i++)
      for (j = 0; j <= i; j++) {
        ker[i][j] = ker[j][i] = kernel(data[i], data[j]);
      }
  }
}
