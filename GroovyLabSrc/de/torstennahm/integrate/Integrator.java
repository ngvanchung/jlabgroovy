/*
 * Created on Oct 20, 2004
 */
package de.torstennahm.integrate;

import java.util.List;

import de.torstennahm.integrate.visualize.Visualizer;

/**
 * The integrator performs integration of the integrand until a specified condition
 * is reached. The condition may for example be that a certain error tolerance
 * has been reached.
 * 
 * The implementing classes use various specific algorithms for numerical
 * integration. Unless stated otherwise, integration is always performed over the unit
 * cube <i>[0,1]^d</i> with the Lesbegue measure, where <i>d</i> is the dimension of the function.
 * <p>
 * Many integrators use quadrature formulas consisting of coordinates and weights.
 * See the package <code>de.torstennahm.integrate.quadratureformula</code>.
 * <p>
 * Implementor's note:
 * Only the abstract method <code>integrate</code> needs to be implemented.
 * All other methods refer to this method using specific stop conditions.
 * In the <code>integrate</code> method, subclasses should perform an
 * open-ended process of integration, in which
 * the result of the integration is continually refined. At regular times, they
 * should check if the integration should be continued. For this, they store the
 * result of the integration so far in an <code>IntegrationResult</code> instance
 * and pass it to the stop condition's <code>stop</code> method.
 * If <code>true</code> is returned, the integration stops
 * and the same <code>IntegrationResult</code> instance for which the
 * stop condition was reached is passed back as final result of the integration. This assures
 * that the <code>IntegrationResult</code> actually fulfills the stop condition.
 * 
 * @see de.torstennahm.integrate.StopConditions
 * 
 * @author Torsten Nahm
 */
public abstract class Integrator<I> {
	/**
	 * Performs numerical integration of the integrand until the given condition
	 * is fulfilled.
	 * 
	 * @param integrand object to be integrated
	 * @param condition stop condition
	 * @param visualizers list of visualizers or <code>null</code> for no visualization
	 * @return result of integration
	 * @throws IntegrationFailedException if an integration error occurs
	 */
	public abstract IntegrationResult integrate(I integrand, StopCondition condition,
			List<Visualizer> visualizers)	throws IntegrationFailedException;

	/**
	 * Convenience method for integration run without visualizers.
	 * Calls <code>integrate(I, StopCondition, null)</code>.
	 * 
	 * @param integrand object to be integrated
	 * @param condition stop condition
	 * @return result of integration
	 * @throws IntegrationFailedException if an integration error occurs
	 */
	public IntegrationResult integrate(I integrand,
			StopCondition condition) throws IntegrationFailedException {
		return integrate(integrand, condition, null);
	}
	
	/**
	 * Convenience method for integration until the specified number of evaluations
	 * has been reached. Calls <code>integrate(integrand, condition, null)</code>
	 * with <code>StopConditions.UntilCallsReached</code> as condition.
	 * 
	 * For some integrators the number of evaluations may exceed
	 * the requested number, but it will never be less than it.
	 * 
	 * @param integrand object to be integrated
	 * @param evals number of function evaluations
	 * @return result of the integration
	 * @throws IntegrationFailedException if an integration error occurs
	 */
	public IntegrationResult integrateByPoints(I integrand,
			long evals) throws IntegrationFailedException {
		return integrate(integrand, new StopConditions.UntilCallsReached(evals));
	}

	/**
	 * Convenience method for integration with visualizers until the specified number of evaluations
	 * has been reached. Calls <code>integrate(integrand, condition, null)</code>
	 * with <code>StopConditions.UntilCallsReached</code> as condition.
	 * 
	 * For some integrators the number of evaluations may exceed
	 * the requested number, but it will never be less than it.
	 * 
	 * @param integrand object to be integrated
	 * @param evals number of evaluations
	 * @return result of the integration
	 * @throws IntegrationFailedException if an integration error occurs
	 */
	public IntegrationResult integrateByPoints(I integrand, long evals,
			List<Visualizer> visualizers) throws IntegrationFailedException {
		return integrate(integrand, new StopConditions.UntilCallsReached(evals));
	}

	/**
	 * Convenience method for integration until the absolute estimated
	 * error is less than than the specified threshold.
	 * Calls <code>integrate(I, StopCondition)</code>
	 * with <code>StopConditions.UntilAbsTol</code> as condition.
	 * Performs numerical integration of the integrand until
	 * 
	 * @param integrand object to be integrated
	 * @param absoluteTolerance absolute tolerance to be reached
	 * @return result of the integration
	 * @throws IntegrationFailedException if an integration error occurs
	 */
	public IntegrationResult integrateAbsTol(I integrand,
			double absoluteTolerance) throws IntegrationFailedException {
		return integrate(integrand, new StopConditions.UntilAbsTol(absoluteTolerance));
	}

	/**
	 * Convenience method for integration with visualizers until the absolute estimated
	 * error is less than than the specified threshold.
	 * Calls <code>integrate(I, StopCondition, List<Visualizer>)</code>
	 * with <code>StopConditions.UntilAbsTol</code> as condition.
	 * Performs numerical integration of the integrand until
	 * 
	 * @param integrand object to be integrated
	 * @param absoluteTolerance absolute tolerance to be reached
	 * @param visualizers list of visualizers or <code>null</code> for no visualization
	 * @return result of the integration
	 * @throws IntegrationFailedException if an integration error occurs
	 */
	public IntegrationResult integrateAbsTol(I integrand, double absoluteTolerance,
			List<Visualizer> visualizers) throws IntegrationFailedException {
		return integrate(integrand, new StopConditions.UntilAbsTol(absoluteTolerance), visualizers);
	}
	
	/**
	 * Convenience method for integration until the absolute estimated
	 * error is less than than the specified threshold.
	 * Calls <code>integrate(I, StopCondition)</code>
	 * with <code>StopConditions.UntilRelTol</code> as condition.
	 * Performs numerical integration of the integrand until
	 * 
	 * @param integrand object to be integrated
	 * @param relativeTolerance absolute tolerance to be reached
	 * @return result of the integration
	 * @throws IntegrationFailedException if an integration error occurs
	 */
	public IntegrationResult integrateRelTol(I integrand,
			double relativeTolerance) throws IntegrationFailedException {
		return integrate(integrand, new StopConditions.UntilRelTol(relativeTolerance));
	}

	/**
	 * Convenience method for integration with visualizers until the relative estimated
	 * error is less than than the specified threshold.
	 * Calls <code>integrate(I, StopCondition)</code>
	 * with <code>StopConditions.UntilRelTol</code> as condition.
	 * Performs numerical integration of the integrand until
	 * 
	 * @param integrand object to be integrated
	 * @param relativeTolerance absolute tolerance to be reached
	 * @param visualizers list of visualizers or <code>null</code> for no visualization
	 * @return result of the integration
	 * @throws IntegrationFailedException if an integration error occurs
	 */
	public IntegrationResult integrateRelTol(I integrand, double relativeTolerance,
			List<Visualizer> visualizers) throws IntegrationFailedException {
		return integrate(integrand, new StopConditions.UntilRelTol(relativeTolerance), visualizers);
	}
	
	/**
	 * Convenience method for integration with visualizers until any one of the conditions
	 * in the list is fulfilled.
	 * Calls <code>integrate(I, StopCondition)</code>
	 * with <code>StopConditions.</code> as condition.
	 * 
	 * @param integrand object to be integrated
	 * @param conditions list of stop conditions
	 * @param visualizers list of visualizers or <code>null</code> for no visualization
	 * @return result of the integration
	 * @throws IntegrationFailedException if an integration error occurs
	 */
	IntegrationResult integrateWithConditions(I integrand, List<StopCondition> conditions,
			List <Visualizer> visualizers) throws IntegrationFailedException {
		return integrate(integrand, new StopConditions.MultiStopCondition(conditions), visualizers);
	}

	/**
	 * Convenience method for integration with visualizers until any one of the conditions
	 * in the list is fulfilled.
	 * Calls <code>integrate(I, StopCondition)</code>
	 * with <code>StopConditions.</code> as condition.
	 * 
	 * @param integrand object to be integrated
	 * @param conditions list of stop conditions
	 * @return result of the integration
	 * @throws IntegrationFailedException if an integration error occurs
	 */
	IntegrationResult integrateWithConditions(I integrand,
			List<StopCondition> conditions) throws IntegrationFailedException {
		return integrate(integrand, new StopConditions.MultiStopCondition(conditions));
	}

}