/*
 * Created on Jul 6, 2003
 */
package de.torstennahm.integrate.quadratureformula;

/**
 * This class generates weights for the trapezoidal quadrature formula on
 * the interval <i>[0,1]</i>.
 * 
 * For <i>n</i> nodes, <i>n>1</i>, the nodes will lie at <i>0</i>,
 * <i>1/(n-1)</i>, <i>2/(n-1)</i>, ..., <i>1</i>.
 * The weights are <i>1/n</i> for all nodes.
 * For <i>n=1</i>, the single <i>(node, weight)</i> pair is <i>(0.5,1)</i>.
 * <p>
 * According to general <code>Generator</code> contract, this class is thread-safe.
 * 
 * @author Torsten Nahm
 */

public class Trapezoidal extends AbstractCachedGenerator {
	public int maxLevel() {
		return 15;
	}
	
	public int maxNodes() {
		return levelToNodes(maxLevel());
	}
	
	@Override
	public QuadratureFormula generateByLevel(int level) {
		if (level < 0 || level > maxLevel()) {
			throw new IllegalArgumentException();
		}
		
		return getByNodes(levelToNodes(level));
	}
	
	private int levelToNodes(int level) {
		if (level == 0) {
			return 1;
		} else {
			return (1 << level) + 1;
		}
	}
	
	@Override
	public QuadratureFormula generateByNodes(int nodesRequested) {
		if (nodesRequested < 1 || nodesRequested > maxNodes()) {
			throw new IndexOutOfBoundsException();
		}
		
		QuadratureFormula qf;
		
		if (nodesRequested == 1) {
			qf = new QuadratureFormula(new double[] {0.5}, new double[] {1.0});
		} else {
			double[] node = new double[nodesRequested], weight = new double[nodesRequested];
			double columnDist = 1.0 / (nodesRequested - 1);
			double columnWeight = 1.0 / (nodesRequested - 1);
			
			for (int i = 0; i < nodesRequested; i++)  {
				node[i] = i * columnDist;
				weight[i] = columnWeight;
			}
			/* Adjust outermost weights. */
			weight[0] *= 0.5;
			weight[nodesRequested - 1] *= 0.5;
			
			qf = new QuadratureFormula(node, weight);
		}
		
		return qf;
	}
	
	@Override
	public String toString() {
		return "Trapezoidal";
	}
}
