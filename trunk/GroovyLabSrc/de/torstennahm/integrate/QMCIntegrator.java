/*
 * Created on Jun 11, 2003
 */
package de.torstennahm.integrate;


import java.util.HashSet;
import java.util.List;
import java.util.Set;

import de.torstennahm.integrate.error.ConvergenceEstimator;
import de.torstennahm.integrate.error.ErrorEstimator;
import de.torstennahm.integrate.visualize.Visualizer;
import de.torstennahm.integrate.visualize.Visualizers;
import de.torstennahm.integrate.visualizerdata.Integrand;
import de.torstennahm.integrate.visualizerdata.NewResult;
import de.torstennahm.integrate.visualizerdata.StartIntegration;
import de.torstennahm.integrate.visualizerdata.StopIntegration;
import de.torstennahm.math.Function;
import de.torstennahm.series.Halton;
import de.torstennahm.series.Series;

/**
 * Performs a Quasi Monte-Carlo integration.
 * 
 * A Halton series is used for the points,
 * of which the first points are discarded. See the documentation of
 * <code>de.torstennahm.math.Halton</code> for details.
 * 
 * @author Torsten Nahm
 */
public class QMCIntegrator extends Integrator<Function> {
	private final PointsGenerator generator;
	
	public QMCIntegrator() {
		this (new HaltonGenerator(Halton.HaltonType.NORMAL));
	}
	
	public QMCIntegrator(PointsGenerator generator) {
		this.generator = generator;
	}
	
	@Override
	public IntegrationResult integrate(Function function, StopCondition condition, List<Visualizer> visualizers) {
		Series<double[]> series = generator.makeSeries(function.inputDimension());
		
		QMCResult result = new QMCResult();
		
		Visualizers.submitToList(visualizers, new Integrand(function));
		Visualizers.submitToList(visualizers, new StartIntegration());
		
		double sum = 0.0;
		int count = 0;
		while (! condition.stop(result)) {
			sum += function.sEvaluate(series.next());
			result.value = sum / ++result.numPoints;
			result.errorEstimator.log(result.numPoints, result.value);
			
			if (++count == 100) {
				Visualizers.submitToList(visualizers, new NewResult(result));
				count = 0;
			}
		}			
		
		Visualizers.submitToList(visualizers, new StopIntegration(result));
		
		return result;
	}
	
	static private class QMCResult implements IntegrationResult {
		private double value = 0.0;
		private long numPoints = 0;
		private ErrorEstimator errorEstimator = new ConvergenceEstimator();
		
		public double value() {
			return value;
		}
		
		public double errorEstimate() {
			return errorEstimator.getEstimate();
		}
		
		public long functionCalls() {
			return numPoints;
		}
		
		public Set<IntegrationInfo> supplementalInfo() {
			return new HashSet<IntegrationInfo>();
		}
	}

	@Override
	public String toString() {
		return "Quasi Monte Carlo Integrator";
	}
}
