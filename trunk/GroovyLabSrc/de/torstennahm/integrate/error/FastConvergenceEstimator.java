/*
 * Created on Mar 29, 2004
 */
package de.torstennahm.integrate.error;

import java.util.LinkedList;
import java.util.ListIterator;

import de.torstennahm.statistics.LinearRegression;

/**
 * Provides an estimate based on the assumption that the integral values
 * shows logarithmic convergence.
 * 
 * By logarithmic convergence, we mean that <i>log(error(n))=a-b*log(n)</i>,
 * where <i>n</i> is the number of function evaluations, and <i>a</i> and
 * <i>b</i> are parameters.
 * 
 * @author Torsten Nahm
 */
public class FastConvergenceEstimator implements ErrorEstimator {
	private double sampleFactor;
	
	private LinkedList<LogEntry> log;
	private long currentPoints, startPoints;
	private double nextPoints;
	private boolean logChanged;
	private double minValue, maxValue;
	
	private double errorStart;
	private double errorSlope;
	
	/**
	 * Create the convergence estimator with a default sample factor of 1.1.
	 * @see #FastConvergenceEstimator(double)
	 */
	public FastConvergenceEstimator() {
		this(1.1);
	}
	
	/**
	 * Create the convergence estimator with the specified sample factor.
	 * The sample factor determines how dense the interpolation points
	 * are to be placed. A factor of <i>f</i> means that when
	 * an interpolation point is set a <i>n</i> function evaluations,
	 * the next should be at <i>f*n</i> evaluations, and so on.
	 * 
	 * @param sampleFactor factor between interpolation points
	 */
	public FastConvergenceEstimator(double sampleFactor) {
		if (sampleFactor <= 1.0) {
			throw new IllegalArgumentException();
		}
		
		this.sampleFactor = sampleFactor;
		
		nextPoints = 100.0;
		log = new LinkedList<LogEntry>();
		minValue = Double.NaN;
		maxValue = Double.NaN;
		startPoints = 0;
		errorStart = Double.NaN;
		errorSlope = Double.NaN;
	}
	
	public void log(long pointsEvaluated, double currentValue) {
		minValue = Math.min(minValue, currentValue);
		maxValue = Math.max(maxValue, currentValue);
		currentPoints = pointsEvaluated;
		
		if (pointsEvaluated >= nextPoints) {
			if (! Double.isNaN(minValue)) {
				LogEntry entry = new LogEntry();
				entry.points = startPoints;
				entry.minValue = minValue;
				entry.maxValue = maxValue;
				
				updateLog(minValue, maxValue);
				
				log.add(entry);
				logChanged = true;
			}
			startPoints = currentPoints;
			nextPoints *= sampleFactor;
			minValue = Double.POSITIVE_INFINITY;
			maxValue = Double.NEGATIVE_INFINITY;
		}
	}
	
	private void updateLog(double minValue, double maxValue) {
		boolean stop = false;
		for (ListIterator<LogEntry> iter = log.listIterator(log.size()); ! stop && iter.hasPrevious(); ) {
			LogEntry entry = iter.previous();
			
			stop = true;
			if (minValue < entry.minValue) {
				entry.minValue = minValue;
				stop = false;
			}
			if (maxValue > entry.maxValue) {
				entry.maxValue = maxValue;
				stop = false;
			}
		}
	}
	
	public double getEstimate() {
		if (logChanged) {
			doRegression();
			logChanged = false;
		}
		return 2 * Math.exp(errorStart + Math.log(currentPoints) * errorSlope);
	}
	
	public double getSlope() {
		return errorSlope;
	}
	
	private void doRegression() {
		errorSlope = Double.NaN;
		errorStart = Double.NaN;
		
		if (currentPoints > 1000) {
			LinearRegression reg = new LinearRegression();
			double logCurrentPoints = Math.log(currentPoints);
			
			for (LogEntry entry : log) {
				double logPoints = Math.log(entry.points);
				if (logPoints >= 0.5 * logCurrentPoints && logPoints <= 0.75 * logCurrentPoints) {
					reg.add(logPoints, Math.log(entry.variation()));
				}
			}
				
			errorSlope = reg.slope();
			errorStart = reg.yIntercept();
		}
	}
	
	static private class LogEntry {
		long points;
		double minValue, maxValue;
		
		double variation() {
			return maxValue - minValue;
		}
	}
}
