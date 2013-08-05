/*
 * Created on Jul 6, 2003
 */
package de.torstennahm.integrate.quadratureformula;

import java.util.HashMap;
import java.util.Map;

/**
 * Provides functionality for caching and retrieving quadrature formulas.
 * <p>
 * Since quadrature formulas are immutable, they can be returned from
 * a cache and need only be calculated once.
 * <p>
 * According to general <code>Generator</code> contract, this class is thread-safe.
 * 
 * @author Torsten Nahm
 */
abstract public class AbstractCachedGenerator implements Generator {
	protected Map<Integer, QuadratureFormula> weightsCache;
	protected Map<Integer, QuadratureFormula> levelCache;
	
	abstract protected QuadratureFormula generateByNodes(int nodesRequested);
	abstract protected QuadratureFormula generateByLevel(int level);
	
	synchronized public QuadratureFormula getByNodes(int nodesRequested) {
		QuadratureFormula w;
		
		if (weightsCache == null) {
			weightsCache = new HashMap<Integer, QuadratureFormula>();
		}
		w = weightsCache.get(new Integer(nodesRequested));
		if (w == null) {
			w = generateByNodes(nodesRequested);
			weightsCache.put(new Integer(nodesRequested), w);
		}
		
		return w;
	}
	
	synchronized public QuadratureFormula getByLevel(int level) {
		QuadratureFormula w;
		
		if (levelCache == null) {
			levelCache = new HashMap<Integer, QuadratureFormula>();
		}
		w = levelCache.get(new Integer(level));
		if (w == null) {
			w = generateByLevel(level);
			levelCache.put(new Integer(level), w);
		}
		
		return w;
	}
}
