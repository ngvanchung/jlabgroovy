/*
 * Created on Jun 11, 2003
 */
package de.torstennahm.integrate;

import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

import de.torstennahm.integrate.visualize.Visualizer;
import de.torstennahm.integrate.visualize.Visualizers;
import de.torstennahm.integrate.visualizerdata.Integrand;
import de.torstennahm.integrate.visualizerdata.NewResult;
import de.torstennahm.integrate.visualizerdata.StartIntegration;
import de.torstennahm.integrate.visualizerdata.StopIntegration;
import de.torstennahm.math.Function;
import de.torstennahm.statistics.Statistics;

/**
 * Performs Monte-Carlo integration, using the java.util.random to generate
 * the points.
 * 
 * The error estimator is based on the assumption that the function's distribution is
 * Gaussian. The error estimate is then multiplied by an arbitrary safety factor of
 * 2 to account for at least some non-Gaussian distributions.
 * 
 * @author Torsten Nahm
 */
public class MCIntegrator extends Integrator<Function> {
	private Random random = new Random();
	
	/*
	 * Given a Gaussian distribution, the chance that the real integral value
	 * lies within the interval given by the methods getResult and getErrorEstimate is
	 * given by the Student T distribution. If after 100 samples the error estimator is to be
	 * be correct with a chance of 0.99, we need the x for which F_T(99)(x) = 0.995,
	 * where F_T(99) is the cumulative distribution function of the T distribution for
	 * n = 99. This value is supplied as T_INTERVAL below. For n > 100, this
	 * value becomes smaller, so that the chance becomes even better than 0.99. It becomes
	 * only slightly smaller, however, so that only little accuracy is sacrified as opposed
	 * to the online calculation of the T values for larger n, which would be rather
	 * expensive.
	 */
	
	/**
	 * Minimum evalutions before error estimate is given.
	 */
	private static final long MIN_EVALUATIONS = 100;
	/**
	 * 99.5%-Quantile of the T(n) distribution for n=99.
	 */
	private static final double T_INTERVAL = 2.625890521;
	
	/**
	 * Arbitraty safety factor that tries to compensate for the fact that the distribution
	 * is in general not Gaussian.
	 */
	private static final double SAFETY_FACTOR = 2.0;

	@Override
	public IntegrationResult integrate(Function function, StopCondition condition, List<Visualizer> visualizers) {
		int dimension = function.inputDimension();
		double[] x = new double[dimension];
		final Statistics statistics = new Statistics();
		
		IntegrationResult result = new IntegrationResult() {
			public double value() {
				return statistics.average();
			}
	
			public double errorEstimate() {
				if (statistics.n() >= MIN_EVALUATIONS) {
					return (statistics.sigma() / Math.sqrt(statistics.n() - 1)) *
						   T_INTERVAL * SAFETY_FACTOR;
				} else {
					return Double.NaN;
				}
			}
	
			public long functionCalls() {
				return (long)statistics.n();
			}
			
			public Set<IntegrationInfo> supplementalInfo() {
				return new HashSet<IntegrationInfo>();
			}
		};
		
		Visualizers.submitToList(visualizers, new Integrand(function));
		Visualizers.submitToList(visualizers, new StartIntegration());
		
		int count = 0;
		while (! condition.stop(result)) {
			for (int j = 0; j < dimension; j++) {
				x[j] = random.nextDouble();
			}
			
			statistics.add(function.sEvaluate(x));
			
			if (++count == 100) {
				Visualizers.submitToList(visualizers, new NewResult(result));
				count = 0;
			}
		}
		
		Visualizers.submitToList(visualizers, new StopIntegration(result));
		
		return result;
	}

	@Override
	public String toString() {
		return "Monte Carlo Integrator";
	}
}
