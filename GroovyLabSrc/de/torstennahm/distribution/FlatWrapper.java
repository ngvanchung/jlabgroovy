/*
 * Created on Nov 30, 2004
 */
package de.torstennahm.distribution;

import de.torstennahm.math.Function;

/**
 * @author Torsten Nahm
 */
public class FlatWrapper extends Function {
	private final Function function;
	private final double center;
	private final double sigma;

	public FlatWrapper(Function function, double center, double variance) {
		this.function = function;
		this.center = center;
		sigma = Math.sqrt(variance);
	}
	
	@Override
	public double sEvaluate(double[] x) {
		double[] xx = x.clone();
		
		double p = center + sigma * x[0];
		xx[0] = p;
		double weight = Math.sqrt(2 * Math.PI) * sigma * Math.exp(p * p / (2 * sigma * sigma));
		
		return weight * function.sEvaluate(xx);
	}

	@Override
	public int inputDimension() {
		return function.inputDimension();
	}
	
	@Override
	public String toString() {
		return "FlatWrapper";
	}
}
