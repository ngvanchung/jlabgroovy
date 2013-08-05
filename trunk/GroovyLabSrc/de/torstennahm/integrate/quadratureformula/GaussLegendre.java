/*
 * Created on Nov 30, 2004
 */
package de.torstennahm.integrate.quadratureformula;

/**
 * This class generates the weights for the Gauss-Legendre quadrature formula.
 * <p> 
 * According to general <code>Generator</code> contract, this class is thread-safe.
 * 
 * @author Torsten Nahm
 */
public class GaussLegendre extends AbstractCachedGenerator {
	private final double EPS = 1e-15;
	
	public int maxLevel() {
		return 8;
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
		if (nodesRequested < 0 || nodesRequested > maxNodes()) {
			throw new IndexOutOfBoundsException();
		}
		
		int n = nodesRequested;
		double[] nodes = new double[n], weights = new double[n];
		
		int m = (n - 1) / 2;
		for (int i = 0; i <= m; i++) {
			double z = Math.cos(Math.PI * (i + 1 - 0.25) / (n + 0.5));
			double pp, jump;
			int iterCount = 0;
			do {
				double p1 = 1.0;
				double p2 = 0.0;
				for (int j = 1; j <= n; j++) {
					double p3 = p2;
					p2 = p1;
					p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j;
				}
				pp = n * (z * p1 - p2) / (z * z - 1.0);
				jump = p1 / pp;
				z -= jump;
				if (++iterCount > 15) {
					throw new RuntimeException("No convergence");
				}
			} while (Math.abs(jump) > EPS);
			
			nodes[i] = 0.5 - 0.5 * z;
			nodes[n - 1 - i] = 0.5 + 0.5 * z;
			weights[i] = 1 / ((1.0 - z * z) * pp * pp);
			weights[n - 1 - i] = weights[i];
		}
		
		return new QuadratureFormula(nodes, weights);
	}
	
	@Override
	public String toString() {
		return "Gauss-Legendre";
	}
}
