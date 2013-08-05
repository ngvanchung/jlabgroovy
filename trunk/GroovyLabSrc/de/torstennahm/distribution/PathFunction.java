/*
 * Created on Aug 29, 2004
 */
package de.torstennahm.distribution;

/**
 * A path function is a weight function for a path integral, with time and place parameters.
 * 
 * For discretized paths, the path cannot be represented as a whole but only the discretized
 * times chosen. This class returns the value of the weighted integral over
 * a (short) path defined only by start and end times, and start and end positions.
 * The class makes a best approximation of the integral with only this information.
 * A general good strategy will be to simply assume a linear path.
 * In many cases, this assumption allows a closed formula for the integral to be found,
 * but the implementation is free in its strategy for giving the approximation.
 * 
 * @author Torsten Nahm
 */
abstract public class PathFunction {
	abstract public double variance();
	abstract public double integrate(double s, double x, double t, double y) ;
}
