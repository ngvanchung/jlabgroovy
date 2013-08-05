/*
 * Created on Mar 17, 2005
 */
package de.torstennahm.integrate;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

import de.torstennahm.integrate.error.ConvergenceEstimator;
import de.torstennahm.integrate.error.ErrorEstimator;
import de.torstennahm.integrate.quadratureformula.Generator;
import de.torstennahm.integrate.sparse.ProductWeightIntegrator;
import de.torstennahm.integrate.visualize.Visualizer;
import de.torstennahm.integrate.visualize.Visualizers;
import de.torstennahm.integrate.visualizerdata.Integrand;
import de.torstennahm.integrate.visualizerdata.StartIntegration;
import de.torstennahm.integrate.visualizerdata.StopIntegration;
import de.torstennahm.math.Function;

/**
 * This type of integrator performs integration by scaling up a set of one-dimensional nodes
 * and weights using the tensor product approach.
 * 
 * Because integration by <code>Integrator</code> subclasses is open-ended, the integrator
 * successively doubles the target number of evaluations and then finds a tensor product grid that
 * approximates this number. This process is repeated until the stop condition is
 * fulfilled.
 * 
 * @author Torsten Nahm
 */
public class ProductIntegrator extends Integrator<Function>  {
	private final ProductWeightIntegrator integrator;
	
	public ProductIntegrator(Generator generator) {
		integrator = new ProductWeightIntegrator(generator);
	}
	
	@Override
	public IntegrationResult integrate(Function integrand, StopCondition condition, List<Visualizer> visualizers) throws IntegrationFailedException {
		int dimension = integrand.inputDimension();
		
		ProductResult result = new ProductResult();
		
		Visualizers.submitToList(visualizers, new Integrand(integrand));
		Visualizers.submitToList(visualizers, new StartIntegration());
		
		while (! condition.stop(result)) {
			int[] nodes = getNodesForEvals(dimension, result.points * 2 + 1);
			result.value = integrator.integrateWithNodes(integrand, nodes);
			result.points += integrator.neededEvaluations(nodes);
		}
		
		Visualizers.submitToList(visualizers, new StopIntegration(result));
		
		return result;
	}
	
	int[] getNodesForEvals(int dimension, long evals) {
		int[] nodes = new int[dimension];
		double nominalNodes = Math.pow(evals, 1.0 / dimension);
		int intNodes = (int) nominalNodes;
		
		double smallEvals = Math.pow(intNodes, dimension);
		int largeNodes = (int) (Math.log(evals/smallEvals) / Math.log((intNodes+1.0)/intNodes) + 0.5);
		for (int i = 0; i < dimension; i++) {
			nodes[i] = i < largeNodes ? intNodes + 1 : intNodes;
		}
		
		return nodes;
	}
	
	static private class ProductResult implements IntegrationResult {
		private double value = 0.0;
		private long points = 0;
		private ErrorEstimator errorEstimator = new ConvergenceEstimator();
		
		public double value() {
			return value;
		}
		
		public double errorEstimate() {
			return errorEstimator.getEstimate();
		}
		
		public long functionCalls() {
			return points;
		}
		
		public Set<IntegrationInfo> supplementalInfo() {
			return new HashSet<IntegrationInfo>();
		}
	}
}
