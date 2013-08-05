/*
 * Created on Jun 3, 2004
 */
package de.torstennahm.integrate.error;

import de.torstennahm.math.MathTN;

/**
 * Estimates the error caused by floating point rounding errors.
 * 
 * The estimator assumes that all calculations are done with
 * <code>double</code> floating point accuracy.
 * 
 * @author Torsten Nahm
 */
public class RoundingErrorEstimator implements ErrorEstimator {
	private double max = 0.0;
	
	public void log(long pointsEvaluated, double currentResult) {
		max = Math.max(max, Math.abs(currentResult));
	}
	
	public double getEstimate() {
		return max * MathTN.FUDGE;
	}
}
