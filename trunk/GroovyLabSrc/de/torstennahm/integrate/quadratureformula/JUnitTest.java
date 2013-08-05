/*
 * Created on Nov 6, 2003
 */
package de.torstennahm.integrate.quadratureformula;

import junit.framework.TestCase;

import de.torstennahm.math.MathTN;

/**
 * @author Torsten Nahm
 */
public class JUnitTest extends TestCase {
	private class TestData {
		Generator g;
		double min, max;
		double weight;
		
		TestData(Generator g, double min, double max, double weight) {
			this.g = g;
			this.min = min;
			this.max = max;
			this.weight = weight;
		}
	}
	
	private TestData[] generatorTest = {
		new TestData(new GaussHermite(), Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, 1.0),
		new TestData(new GaussianWrapper(new Patterson()), Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, 1.0),
		new TestData(new ClenshawCurtis(), 0.0, 1.0, 1.0),
		new TestData(new OpenTrapezoidal(), 0.0, 1.0, 1.0),
		new TestData(new Trapezoidal(), 0.0, 1.0, 1.0),
		new TestData(new Patterson(), 0.0, 1.0, 1.0),
		new TestData(new GaussLegendre(), 0.0, 1.0, 1.0)
	};
	
	/**
	 * Tests the quadrature formula generators.
	 */
	public void testQuadratures() {
		for (int i = 0; i < generatorTest.length; i++) {
			TestData gt = generatorTest[i];
			int max = gt.g.maxNodes();
			if (max == -1 || max >= 1025) {
				max = 1025;
			}
			for (int j = 1; j <= max; j++) {
				QuadratureFormula qf = gt.g.getByNodes(j);
				double totalWeight = 0.0;
				for (int k = 0; k < qf.getSize(); k++) {
					assertTrue(qf.getNode(k) >= gt.min && qf.getNode(k) <= gt.max);
					totalWeight += qf.getWeight(k);
				}
				assertTrue(Math.abs(totalWeight - gt.weight) < MathTN.FUDGE);
			}
		}
	}
	
	/**
	 * Tests the delta wrapper.
	 */
	public void testDeltaQuadrature() {
		for (int i = 0; i < generatorTest.length; i++) {
			TestData gt = generatorTest[i];
			Generator g = new DeltaGenerator(new ClenshawCurtis());
			int max = g.maxLevel();
			if (max == -1 || max > 10) {
				max = 10;
			}
			
			for (int j = 0; j <= max; j++) {
				QuadratureFormula qf = g.getByLevel(j);
				double totalWeight = 0.0;
				for (int k = 0; k < qf.getSize(); k++) {
					assertTrue(qf.getNode(k) >= gt.min && qf.getNode(k) <= gt.max);
					totalWeight += qf.getWeight(k);
				}
				assertTrue(Math.abs(totalWeight - ((j == 0) ? gt.weight : 0.0)) < MathTN.FUDGE);
			}
		}
	}
}
