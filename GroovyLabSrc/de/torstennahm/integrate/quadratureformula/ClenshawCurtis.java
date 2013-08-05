/*
 * Created on Jul 6, 2003
 */
package de.torstennahm.integrate.quadratureformula;

/**
 * This class generates the weights for the Clenshaw-Curtis quadrature formula.
 * <p> 
 * According to general <code>Generator</code> contract, this class is thread-safe.
 * 
 * @author Torsten Nahm
 */

public class ClenshawCurtis extends AbstractCachedGenerator {
	public int maxLevel() {
		return 10;
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
		
		int pairs = nodesRequested;
		if ((pairs & 1) == 0) {		// pairs needs to be odd for Clenshaw-Curtis
			pairs++;
		}
		
		double[] node = new double[pairs], weight = new double[pairs];
		
		double q = 1.0 / (pairs + 1);
		double sum;
		
		for (int i = 1; i <= pairs; i++) {
			node[i - 1] = 0.5 * (1 - Math.cos(Math.PI * i * q));
			sum = 0.0;
			for (int j = 1; j <= pairs; j += 2) {
				sum += Math.sin(j * i * Math.PI / (pairs + 1)) / j;
			}
			weight[i - 1] = 2.0 * q * Math.sin(Math.PI * i * q) * sum;
		}
		
		return new QuadratureFormula(node, weight);
	}
	
	@Override
	public String toString() {
		return "Clenshaw-Curtis";
	}
}
