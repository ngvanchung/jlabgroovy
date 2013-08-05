/*
 * Created on Mar 29, 2004
 */
package de.torstennahm.integrate.error;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import de.torstennahm.statistics.Statistics;

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
public class ConvergenceEstimator implements ErrorEstimator {
	private static final double RESULTSPAN = 100;
	private static final double MINENTRIES = 10;
	
	private double sampleFactor;
	
	private List<LogEntry> log;
	private long currentPoints, startPoints;
	private double nextPoints;
	private boolean logChanged;
	private double estimate = Double.NaN;
	private double minValue, maxValue;
	
	/**
	 * Create the convergence estimator with a default sample factor of 1.1.
	 * @see #ConvergenceEstimator(double)
	 */
	public ConvergenceEstimator() {
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
	public ConvergenceEstimator(double sampleFactor) {
		if (sampleFactor <= 1.0) {
			throw new IllegalArgumentException();
		}
		
		this.sampleFactor = sampleFactor;
		
		nextPoints = 100.0;
		log = new LinkedList<LogEntry>();
		minValue = Double.NaN;
		maxValue = Double.NaN;
		startPoints = 0;
	}
	
	public void log(long pointsEvaluated, double currentValue) {
		minValue = Math.min(minValue, currentValue);
		maxValue = Math.max(maxValue, currentValue);
		currentPoints = pointsEvaluated;
		
		if (pointsEvaluated >= nextPoints) {
			if (! Double.isNaN(minValue)) {
				log.add(0, new LogEntry(startPoints, minValue, maxValue)); 
				logChanged = true;
			}
			startPoints = currentPoints;
			nextPoints *= sampleFactor;
			minValue = Double.POSITIVE_INFINITY;
			maxValue = Double.NEGATIVE_INFINITY;
		}
	}
	
	public double getEstimate() {
		if (logChanged) {
			calcEstimate();
			logChanged = false;
		}
		return estimate;
	}
	
	private void calcEstimate() {
		double[] qB;
		
		qB = calcQB(10.0);
		if (qB != null) {
			qB = calcQB(Math.log(RESULTSPAN + 1) / -Math.log(qB[0]));
		}
		
		if (qB != null) {
			estimate = qB[1] * 2;		// 2 is an arbitrary safety factor
		} else {
			estimate = Double.NaN;
		}
	}
	
	private double[] calcQB(double pointSpan) {
		double[] qB = null;
		double lastResult, lastPoints;
		
		Iterator<LogEntry> iter = log.iterator();
		
		if (iter.hasNext())	{
			LogEntry entry = iter.next();
			lastResult = 0.5 * (entry.minValue + entry.maxValue);
			lastPoints = entry.points;
			
			double totalVar = 0.0;
			boolean spanReached = false;
			int entries = 0;
			int changeEntries = 0;
			while (iter.hasNext() && ! spanReached) {
				entry = iter.next();
				entries++;
				
				double d = Math.max(Math.abs(entry.minValue - lastResult), Math.abs(entry.maxValue - lastResult));
				if (d > totalVar) {
					totalVar = d;
					changeEntries++;
				}
				
				boolean pointSpanReached = lastPoints / entry.points >= pointSpan;
				spanReached = pointSpanReached && (changeEntries >= MINENTRIES);
			}
				
			if (spanReached) {
				String s = null, t = null, u = null;
				long firstPoints = entry.points;
				
				double qMax = 0.0;
				double var = 0;
				Statistics qStat = new Statistics();
				iter = log.iterator();
				iter.next();
				double n = Math.log(lastPoints / firstPoints);
				int entryNum = 0;
				while (entryNum < entries) {
					entry = iter.next();
					
					double d = Math.max(Math.abs(entry.minValue - lastResult), Math.abs(entry.maxValue - lastResult));
					var = Math.max(var, d);
					double m = Math.log(lastPoints / entry.points);
					double qEst = Double.NaN;
					if (var / totalVar < 0.5) {
						qEst = findQ(m , n, var / totalVar);
						if (! Double.isNaN(qEst)) {
							qStat.add(qEst);
							qMax = Math.max(qMax, qEst);
						}
					}
					if (s == null) {
						s = t = u = "[";
					} else {
						s += ", "; t += ", "; u += ", ";
					}
					s += Math.log(entry.points);
					t += qEst;
					u += Math.log(var);
					entryNum++;
				}
				
				double q = (qStat.average() + qStat.sigma() / Math.sqrt(qStat.n() - 1));
				if (0 < q && q < 1) {
					qB = new double[2];
					qB[0] = q;
					qB[1] = totalVar / (Math.pow(q, -n) - 1);
				}
			}
		}
		
		return qB;
	}
	
	private double findQ(double m, double n, double a) {
		double c = m / n;
		if (c < a) {
			return Double.NaN;
		} else {
			double r = 0.0, s = 1.0;
			
			do {
				double t = 0.5 * (r + s);
				double y = (Math.pow(t, -c) - 1.0) / (1.0 / t - 1.0);
				if (y < a) {
					r = t;
				} else {
					s = t;
				}
			} while (Math.abs(r - s) > 0.01);
			
			return Math.pow(0.5 * (r + s), 1.0 / n);
		}
	}
	
	static private class LogEntry {
		long points;
		double minValue, maxValue;
	
		LogEntry(long points, double minValue, double maxValue) {
			this.points = points;
			this.minValue = minValue;
			this.maxValue = maxValue;
		}
	}
}
