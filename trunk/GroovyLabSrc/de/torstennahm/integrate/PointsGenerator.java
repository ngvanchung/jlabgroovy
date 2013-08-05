/*
 * Created on Jan 22, 2005
 */
package de.torstennahm.integrate;

import de.torstennahm.series.Series;

/**
 * Classes implementing this interface produce a <code>Series</code> of points
 *  These points could for example be generated randomly
 * (Monte Carlo) or be from a Halton series (for Quasi Monte Carlo). All points
 * must lie in the unit cube <i>[0,1]^d</i>, where <i>d</i> is the given dimension.
 * 
 * @author Torsten Nahm
 */
public interface PointsGenerator {
	/**
	 * Creates a <code>Series</code> object, which supplies the points.
	 * 
	 * @param dimension dimension of the unit cube from which the points are selected
	 * @return <code>Series</code> iterator
	 * @see de.torstennahm.series.Series
	 */
	Series<double[]> makeSeries(int dimension);
}
