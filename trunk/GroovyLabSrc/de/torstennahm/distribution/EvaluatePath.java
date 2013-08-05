/*
 * Created on Jul 20, 2003
 */
package de.torstennahm.distribution;

import de.torstennahm.math.Function;


/**
 * Integrates an integrand function over a path.
 * 
 * The given function <i>f</i> must have dimension 2, with the first argument
 * being the time, and the second the path coordinate at the time. For this function, the
 * integral <i>f(t,x) dt</i> over the interval <i>[tStart, tEnd]</i> is evaluated
 * numerically with the trapezoid rule.
 * <p>
 * According to general <code>Function</code> contract, this class is thread-safe.
 * 
 * @author Torsten Nahm
 */
public class EvaluatePath extends Function {
	protected final int steps;
	protected final double tStart, tEnd;
	protected final Function integrand;
	
	/**
	 * Constructs the path evaluator, specifying the start time, end time and steps
	 * of the path. Note that the case <i>tStart > tEnd</i> is explicitly supported.
	 * The integral is not oriented, that is the integral over the constant 1
	 * will the <i>abs(tStart - tEnd)</i> in both cases.
	 * 
	 * @param integrand Function with dimension 2
	 * @param tStart start time of the path
	 * @param tEnd end time of the path
	 * @param steps number of steps for the path, must be >= 1
	 * @throws IllegalArgumentException if integrand does not have dimension 2
	 * @throws IllegalArgumentException if steps < 1
	 * @throws NullPointerException if integrand is null
	 */
	public EvaluatePath(Function integrand, double tStart, double tEnd, int steps) {
		if (integrand.inputDimension() != 2) {
			throw new IllegalArgumentException("Integrand must have dimension 2");
		}
		
		this.integrand = integrand;
		this.tStart = tStart;
		this.tEnd = tEnd;
		this.steps = steps;
	}
	
	@Override
	public int inputDimension() {
		return steps + 1;
	}
	
	/**
	 * Evaluates the integral of the integrand over the given path <code>x</code>.
	 * <code>x</code> represents the path as follows: <code>x[0]</code> is
	 * the position of the path at the time <i>tStart</i>, <code>x[1]</code>
	 * the position at <i>tStart + (tEnd - tStart) / steps</i>, ...,
	 * <code>x[steps]</code> the position at <i>tEnd</i>.
	 * Integration is performed using the trapezoid rule.
	 * 
	 * @param x the positions of the path at equidistant time steps
	 */
	@Override
	public double sEvaluate(double[] x) {
		checkArgument(x);
		
		double sum;
		double[] tx = new double[2];
		double dT = (tEnd - tStart) / steps;
		
		/* Calculate the leftmost and rightmost value seperately, so the inner loop
		 * is not slowed down by having to use different weights
		 */
		tx[0] = tStart;
		tx[1] = x[0];
		sum = integrand.sEvaluate(tx);
		tx[0] = tEnd;
		tx[1] = x[steps];
		sum += integrand.sEvaluate(tx);
		sum *= 0.5;
		
		tx[0] = tStart;
		for (int i = 1; i < steps; i++) {
			tx[0] += dT;
			tx[1] = x[i];
			sum += integrand.sEvaluate(tx);
		}
		
		return sum * Math.abs(tEnd - tStart) / steps;
	}
	
	@Override
	public String toString() {
		return "Path integral for " + integrand;
	}
}
