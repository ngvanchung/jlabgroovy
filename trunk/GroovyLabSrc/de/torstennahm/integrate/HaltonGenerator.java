/*
 * Created on Jan 27, 2005
 */
package de.torstennahm.integrate;

import de.torstennahm.series.Halton;
import de.torstennahm.series.Series;
import de.torstennahm.series.Halton.HaltonType;

/**
 * Generates a Halton series.
 * 
 * @author Torsten Nahm
 */
public class HaltonGenerator implements PointsGenerator {
	private final HaltonType type;
	
	/**
	 * Construct the generator.
	 * 
	 * @param type type of Halton series as defined in @see de.torstennahm.series.Halton
	 */
	public HaltonGenerator(Halton.HaltonType type) {
		this.type = type;
	}
	
	public Series<double[]> makeSeries(int dimension) {
		return new Halton(dimension, type);
	}
}
