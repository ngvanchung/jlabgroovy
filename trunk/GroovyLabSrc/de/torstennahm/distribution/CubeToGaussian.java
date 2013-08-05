/*
 * Created on Feb 4, 2004
 */
package de.torstennahm.distribution;

import de.torstennahm.math.MathTN;
import de.torstennahm.math.VectorFunction;

/**
 * Creates Gaussian output from input on the multi-dimensional unit cube.
 * 
 * Specifically, the class converts input with the <i>n</i>-dimensional
 * uniform distribution on the unit cube to output with 
 * an <i>n</i>-dimensional Gaussian distibution.
 * 
 * @author Torsten Nahm
 */
public class CubeToGaussian extends VectorFunction {
	protected final int dimension;
	
	public CubeToGaussian(int dimension) {
		this.dimension = dimension;
	}
	
	@Override
	public int inputDimension() {
		return dimension;
	}
	
	@Override
	public int outputDimension() {
		return dimension;
	}

	@Override
	public double[] evaluate(double[] x) {
		checkArgument(x);
		
		double[] gx = new double[x.length];
		
		for (int i = 0; i < gx.length; i++) {
			gx[i] = MathTN.inverseGaussian(x[i]);
		}
		
		return gx;
	}
	
	@Override
	public String toString() {
		return "Cube to Gauss";
	}
}
