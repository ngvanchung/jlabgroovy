/*
 Created on Jul 20, 2003
 */
package de.torstennahm.distribution;

import de.torstennahm.math.Function;


/**
 * Numerically evaluates the Feynman-Kac formula for an integrand and border
 * function over a path.
 * 
 * It numerically evaluates <i>b(x(t)) * exp(integral_r=0..t v(t-r,x(r))) dr)</i>,
 * where <i>b</i> is the border function, <i>v</i> is the integrand function,
 * and <i>x</i> is the path.
 * <p>
 * According to general <code>Function</code> contract, this class is thread-safe.
 * 
 * @author Torsten Nahm
 */
public class FeynmanKac extends EvaluatePath {
	protected final Function border;
	
	/**
	 * Constructs the Feynman-Kac object.
	 * 
	 * @param integrand Function to be integrated over the path. It must have dimension 2,
	 * with the first argument being the time, and the second the path position.
	 * @param border Border function. Its argument is the end position of the path.
	 * @param timeSpan time span of the path
	 * @param steps number of steps of the path
	 * @throws NullPointerException if integrand or border is null
	 */
	public FeynmanKac(final Function integrand, Function border, final double timeSpan, int steps) {
		super(integrand, 0.0, timeSpan, steps);
		
		if (border == null) {
			throw new NullPointerException("Border may not be null");
		}
		
		this.border = border;
	}
	
	/**
	 * Evaluates the integral of the integrand over the given path <code>x</code>.
	 * <code>x</code> represents the path as follows: <code>x[0]</code> is
	 * the position of the path at the time <i>0</i>, <code>x[1]</code>
	 * the position at <i>timeSpan / steps</i>, ...,
	 * <code>x[steps]</code> the position at <i>timeSpan</i>.
	 * 
	 * @param x the positions of the path at equidistant time steps
	 */
	@Override
	public double sEvaluate(double[] x) {
		checkArgument(x);
		return border.sEvaluate(new double[] {x[0]}) * Math.exp(super.sEvaluate(x));
	}
	
	@Override
	public String toString() {
		return "Feynman-Kac for the integrand " + integrand + " and the border " + border;
	}
}
