/*
 * Created on Jul 6, 2003
 */
package de.torstennahm.integrate.quadratureformula;


/**
 * This class generates weights for the modified trapezoidal quadrature formula for
 * the open interval ]0,1[.
 * 
 * For <i>n</i> nodes, <i>n>1</i>, the nodes will lie at
 * <i>1/(n+1)</i>, <i>2/(n+1)</i>, ..., <i>n/(n+1)</i>.
 * The weights are <i>1/(n+1)</i> for all internal nodes and
 * <i>1.5/(n+1)</i> for the first and the last node.
 * For <i>n=1</i>, the single <i>(node, weight)</i> pair is <i>(0.5,1)</i>.
 * <p>
 * According to general <code>Generator</code> contract, this class is thread-safe.
 * 
 * @author Torsten Nahm
 */

public class OpenTrapezoidal extends AbstractCachedGenerator {
	public int maxLevel() {
		return 15;
	}
	
	public int maxNodes() {
		return levelToNodes(maxLevel());
	}
	
	private int levelToNodes(int level) {
		return (1 << (level + 1)) - 1;
	}
	
	@Override
	public QuadratureFormula generateByLevel(int levelRequested) {
		if (levelRequested < 0 || levelRequested > maxLevel()) {
			throw new IllegalArgumentException();
		}
		
		return getByNodes(levelToNodes(levelRequested));
	}
	
	@Override
	public QuadratureFormula generateByNodes(int nodesRequested) {
		if (nodesRequested < 1 || nodesRequested > maxNodes()) {
			throw new IllegalArgumentException();
		}
		
		QuadratureFormula qf;
		
		if (nodesRequested == 1) {
			qf = new QuadratureFormula(new double[] {0.5}, new double[] {1.0});
		} else {
			double[] node = new double[nodesRequested], weight = new double[nodesRequested];
			double columnWeight = 1.0 / (nodesRequested + 1);
			double columnDist = 1.0 / (nodesRequested + 1);
			
			for (int i = 0; i < nodesRequested; i++)  {
				node[i] = (i + 1) * columnDist;
				weight[i] = columnWeight;
			}
			/* Adjust outermost weights for open integration. */
			weight[0] *= 1.5;
			weight[nodesRequested - 1] *= 1.5;
			
			qf = new QuadratureFormula(node, weight);
		}
		
		return qf;
	}
	
	@Override
	public String toString() {
		return "Open Trapezoidal";
	}
}
