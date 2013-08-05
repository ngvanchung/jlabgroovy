/*
 * Created on Jul 6, 2003
 */
package de.torstennahm.integrate.quadratureformula;

import de.torstennahm.math.MathTN;

/**
 * This generator takes an underlying quadrature formula generator for the uniform
 * distribution on the interval <i>[0,1]</i>, producing a generator for the Gaussian
 * normal distribution.
 * 
 * Using this new quadrature formula generator
 * corresponds to integrating with the normal distribution.
 * <p>
 * The wrapper works by transforming the nodes of the quadrature formula
 * with the inverse of the normal Gaussian cumulative distribution function.
 * <p>
 * According to general <code>Generator</code> contract, this class is thread-safe.
 *  
 * @author Torsten Nahm
 */
public class GaussianWrapper extends AbstractCachedGenerator {
	final private Generator generator;
	
	/**
	 * Creates the <code>InverseGaussianGenerator</code>.
	 * 
	 * @param generator The underlying generator
	 */
	public GaussianWrapper(Generator generator) {
		this.generator = generator;
	}
	
	public int maxLevel() {
		return generator.maxLevel();
	}
	
	public int maxNodes() {
		return generator.maxNodes();
	}
	
	@Override
	public QuadratureFormula generateByNodes(int nodesRequested) {
		return transform(generator.getByNodes(nodesRequested));
	}
	
	@Override
	public QuadratureFormula generateByLevel(int level) {
		return transform(generator.getByLevel(level));
	}
	
	private QuadratureFormula transform(QuadratureFormula oldQF) {
		int size = oldQF.getSize();
		double[] node = new double[size], weight = new double[size];
		for (int i = 0; i < size; i++) {
			node[i] = MathTN.inverseGaussian(oldQF.getNode(i));
			weight[i] = oldQF.getWeight(i);
		}
		
		return new QuadratureFormula(node, weight); 
	}
	
	@Override
	public String toString() {
		return "Gaussian " + generator.toString();
	}
}
