/*
 * Created on Jan 22, 2007
 */
package de.torstennahm.distribution;

/**
 * Models a 1-dimensional probability kernel of a Markov process.
 * 
 * @author Torsten Nahm
 */
public interface Kernel {
	/**
	 * Returns a new random value for the process based on the distribution
	 * of the kernel.
	 * 
	 * @param time time span
	 * @param start starting position
	 * @param random a random number with uniform distribution in [0,1)
	 * @return new position for the particle after the given time
	 */
	double next(double time, double start, double random);
}
