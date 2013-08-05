/*
 * Created on Jan 22, 2007
 */
package de.torstennahm.distribution;

import de.torstennahm.math.MathTN;

/**
 * Implements a kernel for the normal Gaussian distribution with variance 1 (<i>N(0,1)</i>).
 * 
 * @author Torsten Nahm
 */
public class GaussianKernel implements Kernel {
	public double next(double time, double start, double random) {
		return start + Math.sqrt(time) * MathTN.inverseGaussian(random);
	}
}
