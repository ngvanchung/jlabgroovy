/*
 * Created on May 30, 2004
 */
package de.torstennahm.integrate.error;

/**
 * Provides an estimate of the current error of the integral.
 * 
 * The estimate is based on a log of previous integration values.
 * 
 * @author Torsten Nahm
 */
public interface ErrorEstimator {
	/**
	 * Logs the current integration value.
	 * 
	 * @param pointsEvaluated number of function evaluations for this integral value
	 * @param currentValue current integral value
	 */
	public abstract void log(long pointsEvaluated, double currentValue);
	
	/**
	 * Returns the estimate of the integration error based on the log.
	 * 
	 * @return estimate of the integration error
	 */
	public abstract double getEstimate();
}