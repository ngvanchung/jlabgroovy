/*
 * Created on Oct 20, 2004
 */
package de.torstennahm.integrate.sparse;

import java.util.List;

import de.torstennahm.integrate.Integrator;
import de.torstennahm.integrate.IntegrationFailedException;
import de.torstennahm.integrate.IntegrationResult;
import de.torstennahm.integrate.StopCondition;
import de.torstennahm.integrate.quadratureformula.Generator;
import de.torstennahm.integrate.quadratureformula.Patterson;
import de.torstennahm.integrate.sparse.evaluateindex.DeltaWeightEvaluator;
import de.torstennahm.integrate.sparse.evaluateindex.Evaluator;
import de.torstennahm.integrate.visualize.Visualizer;
import de.torstennahm.math.Function;

/**
 * A convenience class that uses generally acceptable defaults for integration.
 * It performs adaptive sparse grid integration for the uniform measure
 * on <i>[0,1]^d</i> using <code>Patterson</code> weights and the
 * <code>EvaluateIntegrator</code> strategy.
 * 
 * If you have additional information about the functions, you should consider using a more
 * specific integrator suited to the particular function class in question. 
 * 
 * @see de.torstennahm.integrate.sparse.EvaluateIntegrator
 * @see de.torstennahm.integrate.quadratureformula.Patterson
 * 
 * @author Torsten Nahm
 */
public class DefaultSparseIntegrator extends Integrator<Function> {
	private final Integrator<Evaluator> integrator = new EvaluateIntegrator();
	private final Generator generator;
	
	/**
	 * Simple, optionless constructor.
	 */
	public DefaultSparseIntegrator() {
		this (new Patterson());
	}
	
	/**
	 * Constructor that allows specification of the integration formula.
	 * 
	 * @param generator generator for the weights and nodes
	 */
	public DefaultSparseIntegrator(Generator generator) {
		this.generator = generator;
	}
	
	@Override
	public IntegrationResult integrate(Function function, StopCondition condition, List<Visualizer> visualizers) throws IntegrationFailedException {
		return integrator.integrate(new DeltaWeightEvaluator(function, generator), condition, visualizers);
	}
}
