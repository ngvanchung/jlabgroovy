/*
 * Created on Mar 29, 2004
 */
package de.torstennahm.integrate.error;

/**
 * Provides a very simple estimate of the integration error.
 * 
 * The estimate given is the absolute difference between the integral
 * value at <i>n/2</i> and <i>n</i> function evaluations, where
 * <i>n</i> increases with new log entries.
 * 
 * @author Torsten Nahm
 */
public class SimpleEstimator implements ErrorEstimator {
	private double lastResult = Double.NaN;
	private double estimate = Double.NaN;
	private long nextStop = 100;
	
	public void log(long pointsEvaluated, double currentResult) {
		if (pointsEvaluated >= nextStop) {
			estimate = Math.abs(lastResult - currentResult);
			lastResult = currentResult;
			nextStop *= 2;
		}
	}
	
	public double getEstimate() {
		return estimate;
	}
}
