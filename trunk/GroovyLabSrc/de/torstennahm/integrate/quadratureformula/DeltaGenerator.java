/*
 * Created on Jul 6, 2003
 */
package de.torstennahm.integrate.quadratureformula;

import java.util.*;

import de.torstennahm.math.MathTN;

/**
 * This generator takes an underlying quadrature formula generator,
 * producing a corresponding generator for delta integration.
 * 
 * Let <code>wrapped</code> be the underlying generator.
 * Then, for or a given level <code>l</code>, <code>this.getByLevel(l)</code> will return
 * a quadrature formula representing the difference between <code>wrapped.getByLevel(l)</code>
 * minus <code>wrapped.getByLevel(l-1)</code>.
 * In the special case of <code>l</code> equals 0, <code>wrapped.getByLevel(0)</code> is
 * returned.
 * <p>
 * The <code>DeltaWrapper</code> generator supports only the <code>getByLevel</code> method,
 * as the <code>getByNodes</code> method does not make sense for delta integration.
 * <p>
 * According to general <code>Generator</code> contract, this class is thread-safe.
 * 
 * @author Torsten Nahm
 */
public class DeltaGenerator extends AbstractCachedGenerator {
	final private Generator generator;
	
	/**
	 * Creates the <code>DeltaGenerator</code>.
	 * 
	 * @param generator generator for which the delta quadrature formula are to be calculated
	 */
	public DeltaGenerator(Generator generator) {
		this.generator = generator;
	}
	
	public int maxLevel() {
		return generator.maxLevel();
	}
	
	/**
	 * This method is not supported, as delta integration only makes
	 * sense between successive levels of quadrature formulas.
	 * 
	 * @throws UnsupportedOperationException
	 */
	@Override
	public QuadratureFormula generateByNodes(int nodesRequested) {
		throw new UnsupportedOperationException();
	}
	
	/**
	 * This method is not supported, as delta integration only makes
	 * sense between successive levels of quadrature formulas.
	 * 
	 * @throws UnsupportedOperationException
	 */
	public int maxNodes() {
		throw new UnsupportedOperationException();
	}
	
	/**
	 * Returns the quadrature formula representing the difference between
	 * quadrature formula produced by the underlying generator for the specified level
	 * and the specified level minus 1. If <code>level</code> is <code>0</code>,
	 * the unmodified weights for level 0 are returned.
	 * 
	 * @return delta quadrature formula for the specified level
	 */
	@Override
	public QuadratureFormula generateByLevel(int level) {
		QuadratureFormula w0, w1;
		QuadratureFormula w;
		
		w1 = generator.getByLevel(level);
		if (w1 == null) {
			w = null;
		} else if (level == 0) {
			w = w1;
		} else {
			w0 = generator.getByLevel(level - 1);
			
			// AUTOBOXING
			List<Double> nodeList = new LinkedList<Double>(), weightList = new LinkedList<Double>();
			int i0 = 0, i1 = 0;
			int s0 = w0.getSize(), s1 = w1.getSize();
			
			/* Consolidate the weights, adding weights with the same coordinates together */
			
			while (i0 < s0 || i1 < s1) {
				if (i0 < s0 && i1 < s1 && Math.abs(w0.getNode(i0) - w1.getNode(i1)) < MathTN.FUDGE) {
					nodeList.add(new Double(w0.getNode(i0)));
					weightList.add(new Double(-w0.getWeight(i0) + w1.getWeight(i1)));
					i0++;
					i1++;
				} else if (i1 == s1 || (i0 < s0 && w0.getNode(i0) < w1.getNode(i1))) {
					nodeList.add(new Double(w0.getNode(i0)));
					weightList.add(new Double(-w0.getWeight(i0)));
					i0++;
				} else {
					nodeList.add(new Double(w1.getNode(i1)));
					weightList.add(new Double(w1.getWeight(i1)));
					i1++;
				}
			}
			
			w = new QuadratureFormula(nodeList, weightList);
		}
		
		return w;
	}
	
	@Override
	public String toString() {
		return "Delta " + generator.toString();
	}
}
