/*
 * Created on Mar 29, 2004
 */
package de.torstennahm.integrate.error;

import java.util.LinkedList;
import java.util.List;

import de.torstennahm.math.MathTN;

/**
 * Models the integration as a random walk with decreasing variation per step.
 * 
 * The maximum likelihood estimator for the error under this model is used.
 * 
 * @author Torsten Nahm
 */
public class RandomWalkEstimator implements ErrorEstimator {
	private int samples;
	private double sampleFactor;
	private List<Double> log;
	private double nextPoints;
	private boolean logChanged;
	private double estimate;
	
	public RandomWalkEstimator() {
		this(15, 1.1);
	}
	
	public RandomWalkEstimator(int samples, double sampleFactor) {
		if (samples <= 2 || sampleFactor <= 1.0) {
			throw new IllegalArgumentException();
		}
		
		this.samples = samples;
		this.sampleFactor = sampleFactor;
		
		nextPoints = 100.0;
		log = new LinkedList<Double>();
	}
	
	public void log(long pointsEvaluated, double currentResult) {
		if (pointsEvaluated >= nextPoints) {
			log.add(0, new Double(currentResult)); 
			if (log.size() > samples) {
				log.remove(samples);
			}
			logChanged = true;
			nextPoints *= sampleFactor;
		}
	}
	
	public double getEstimate() {
		if (log.size() < samples) {
			return Double.NaN;
		} else {
			if (logChanged) {
				calcEstimate();
				logChanged = false;
			}
			return estimate;
		}
	}
	
	private void calcEstimate() {
		double[] polynomial = new double[samples - 1];
		double[] d = new double[samples];
		
		for (int i = 0; i < samples; i++) {
			d[i] = log.get(i).doubleValue();
		}
		
		for (int i = 1; i < samples; i++) {
			double c = -samples + 2 * i;
			polynomial[i - 1] = c * (d[i] - d[i-1]) * (d[i] - d[i-1]);
		}
		
		double q = solveNewton(polynomial);
		
		double sum = 0.0;
		double qPow = 1.0;
		for (int i = 1; i < samples; i++) {
			sum += (d[i] - d[i-1]) * (d[i] - d[i-1]) * qPow;
			qPow *= q;
		}
		
		double b = sum / (samples - 1) * q / (1 - q);
		
		estimate = MathTN.inverseGaussian(0.99) * Math.sqrt(b);
	}
	
	private double solveNewton(double[] p) {
		double x = 0.9, xx;
		int count = 20;
		
		double[] pDiff = new double[p.length - 1];
		for (int i = 0; i < pDiff.length; i++) {
			pDiff[i] = p[i + 1] * (i + 1);
		}
		
		do {
			xx = x;
			x -= evalPolynomial(p, x) / evalPolynomial(pDiff, x);
		} while(Math.abs(x - xx) > MathTN.FUDGE && count-- > 0);
		
		return (count == 0) ? Double.NaN : x;
	}
	
	private double evalPolynomial(double[] p, double x) {
		double r = 0.0;
		double pot = 1.0;
		
		for (int i = 0; i < p.length; i++) {
			r += p[i] * pot;
			pot *= x;
		}
		
		return r;
	}
}
