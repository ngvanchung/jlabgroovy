/*
 * Created on Jul 6, 2003
 */
package de.torstennahm.integrate.quadratureformula;


/**
 * This class generates quadrature formulas for given levels of integration. Subclasses
 * each implement specific rules, such as the trapezoidal rule or the Clenshaw-Curtis rule.
 * <p>
 * Implementor's note: The generation of the quadrature formulas should be quick.
 * It is therefore recommended to use a cache for the quadrature formulas, so
 * they can be returned immediately. This may be done by subclassing
 * <code>AbstractGeneratorCache</code>, which provides a framework for caching.
 * <p>
 * All subclasses are required to be thread-safe.
 * 
 * @author Torsten Nahm
 */
public interface Generator {
	/**
	 * Returns a set of node-and-weight pairs in an <code>QuadratureFormula</code>
	 * object. The number of nodes returned should be the same as those requested.
	 * If this is not possible due to the special structure of the generator, the number
	 * of nodes returned must be larger than that requested.
	 * However, when 1 node is requested exactly 1 pair must be returned.
	 * <p>
	 * The number of nodes requested must be between 1 and <code>getMaxNodes()</code>.
	 * 
	 * @param nodesRequested number of nodes for the quadrature formula
	 * @return quadrature formula
	 * @throws IllegalArgumentException if the number of nodes requested is not in the valid range
	 */
	QuadratureFormula getByNodes(int nodesRequested);
	
	/**
	 * This routine is similar to <code>getByLevel</code>, but here the argument
	 * is only an abstract measure of the number of node-and-weight-pairs returned.
	 * In general, the number of weights should rise exponentially with the level.
	 * However, the only conditions imposed it that the number of pairs should rise
	 * strictly monotonously with the level. <code>level</code> must be between
	 * 0 and <code>getMaxLevel()</code>.
	 * <p>
	 * In the special case where <code>level</code> is requested as 0, a quadrature formula
	 * with exactly 1 pair must be returned.
	 * 
	 * @param level level of the quadrature formula
	 * @return quadrature formula
	 * @throws IllegalArgumentException if the level is not in the valid range
	 */
	QuadratureFormula getByLevel(int level);
	
	/**
	 * Returns the maximum argument supported by <code>getByNodes</code>.
	 * 
	 * @return maximum number of nodes, or -1 for no maximum
	 */
	int maxNodes();
	
	/**
	 * Returns the maximum argument supported by <code>getByLevel</code>.
	 * 
	 * @return maximum level, or -1 for no maximum
	 */
	int maxLevel();
}
