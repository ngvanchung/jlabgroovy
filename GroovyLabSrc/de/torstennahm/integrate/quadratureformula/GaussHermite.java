/*
 * Created on Jul 6, 2003
 */
package de.torstennahm.integrate.quadratureformula;

import de.torstennahm.util.Util;
import flanagan.roots.RealRoot;
import flanagan.roots.RealRootDerivFunction;

/**
 * Generates the Gauss-Hermite quadrature formula.
 * 
 * The quadrature formulas approximate a Gaussian distribution on the real
 * numbers with the specified sigma.
 * <p>
 * The quadrature formula is calculated using an iterative approximation.
 * <p>
 * According to general <code>Generator</code> contract, this class is thread-safe.
 * 
 * @author Torsten Nahm
 */

public class GaussHermite extends AbstractCachedGenerator {
	private final double sigma;
	private final double factor;
	static private final double PiRoot = Math.pow(Math.PI, -1.0/2);
	static private final double PiRoot4 = Math.pow(Math.PI, -1.0/4);
	
	/**
	 * Constructs the generator for the Gaussian distribution with sigma
	 * equal to sqrt(2).
	 */
	public GaussHermite() {
		this(Math.sqrt(2.0));
	}
	
	/**
	 * Constructs the generator for the Gaussian distribution with the 
	 * specified sigma.
	 * 
	 * @param sigma sigma of the requested Gaussian distribution
	 */
	public GaussHermite(double sigma) {
		this.sigma = sigma;
		factor = Math.sqrt(2.0) / sigma;
	}
	
	public int maxLevel() {
		int level = 0;
		
		while (levelToNodes(level + 1) <= maxNodes()) {
			level++;
		}
		
		return level;
	}
	
	public int maxNodes() {
		return 400;		// More nodes lead to floating point overflows/underflows
	}
	
	private int levelToNodes(int level) {
		return (1 << (level + 1)) - 1;
	}
	
	@Override
	protected QuadratureFormula generateByLevel(int levelRequested) {
		if (levelRequested < 0 || levelRequested > maxLevel()) {
			throw new IllegalArgumentException("Level not in valid range");
		}
		
		return getByNodes(levelToNodes(levelRequested));
	}
	
	@Override
	protected QuadratureFormula generateByNodes(int nodesRequested) {
		if (nodesRequested < 0 || nodesRequested > maxNodes()) {
			throw new IllegalArgumentException("Nodes requested not in valid range");
		}
		
		double[] nodes, weights;
		
		nodes = new double[nodesRequested];
		weights = new double[nodesRequested];
		
		if (nodesRequested == 1) {
			nodes[0] = 0.0;
			weights[0] = 1.0;
		} else {
			int n = nodesRequested;
			int m = (n + 1) / 2;
			double tol = 1e-14;
			
			HermiteFunction hFunction = new HermiteFunction(n);
			
			double z = 0.0;
			
			for (int i = 0; i < m; i++) {
				if (i == 0) {
					z =	Math.sqrt((2 * n + 1)) - 1.85575 * Math.pow((2 * n + 1), -1.0/6);
					z = RealRoot.newtonRaphson(hFunction, z, tol);
				} else if (i == 1) {
					double lowZ = z - 1.5 * Math.pow(n, 0.426) / z;
					z = RealRoot.bisectNewtonRaphson(hFunction, lowZ, z - tol, tol);
				} else {
					double lowZ = z + (z - nodes[i - 2]);
					z = RealRoot.bisectNewtonRaphson(hFunction, lowZ, z - tol, tol);
				}
				
				double pp = hFunction.function(z)[1];
				
				nodes[i] = z;
				nodes[n - 1 - i] = -nodes[i];
				weights[n - 1 - i] = weights[i] = PiRoot * 2.0 / (pp * pp);
			}
			
			for (int i = 0; i < n; i++) {
				nodes[i] *= factor;
			}
			nodes = Util.reverseArray(nodes);
		}
		
		return new QuadratureFormula(nodes, weights);
	}
	
	@Override
	public String toString() {
		return "Gauss-Hermite (sigma=" + sigma + ")";
	}
	
	private static class HermiteFunction implements RealRootDerivFunction {
		private final int n;
		
		HermiteFunction(int n) {
			this.n = n;
		}
		
		public double[] function(double z) {
			double p1 = PiRoot4;
			double p2 = 0.0;
			for (int j = 0; j < n; j++) {
				double p3 = p2;
				p2 = p1;
				p1 = z * Math.sqrt(2.0 / (j + 1)) * p2 - Math.sqrt((double)j / (j + 1)) * p3;
			}
			
			return new double[] { p1, Math.sqrt(2.0 * n) * p2 };
		}
	}
}
