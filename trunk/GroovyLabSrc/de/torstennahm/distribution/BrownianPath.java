/*
 * Created on Oct 12, 2004
 */
package de.torstennahm.distribution;

import de.torstennahm.math.VectorFunction;
import de.torstennahm.util.Util;

/**
 * Creates Brownian paths from Gassian input using the the Brownian bridge
 * algorithm.
 * This class is built on top of the <code>BrownianBridge</code> and <code>BrownianPath</code>
 * classes, and offers additional functionality.
 * In particular, it allows open, closed and trace paths.
 * 
 * @author Torsten Nahm
 * 
 * This class is thread-safe.
 */
public class BrownianPath extends VectorFunction {
	static public final int BRIDGE = 0, WALK = 1;
	static private final int FIXED_AT_START = 0, FIXED_AT_END = 1, FIXED = 2, TRACE = 3;
	
	private int steps = 2;
	private int generatorType = BRIDGE;
	private int type = FIXED_AT_START;
	private double variance = 1.0;
	private double start = 0.0, end;
	
	private int endDimensions;
	private VectorFunction pathFunction;
	
	/**
	 * Set to number of steps to the given value. The output dimension of the function
	 * is <code>steps+1</code>. The default number of steps is 2.
	 * 
	 * @param steps number of steps
	 */
	public synchronized void setSteps(int steps) {
		this.steps = steps;
		pathFunction = null;
	}
	
	/**
	 * Sets the method by which the Brownian path is created, either
	 * <code>BrownianPath.BRIDGE</code> for a Brownian bridge or <code>BrownianPath.WALK</code>
	 * for a random walk. Note that closed paths are only possible for a Brownian bridge.
	 * The default is <code>BrownianPath.BRIDGE</code>.
	 * 
	 * @param type 
	 */
	public synchronized void setType(int type) {
		generatorType = type;
		pathFunction = null;
	}
	
	/**
	 * 
	 * @param start
	 * @param variance
	 */
	public synchronized void makeFixedAtStart(double start, double variance) {
		type = FIXED_AT_START;
		this.start = start;
		this.variance = variance;
		endDimensions = 1;
		pathFunction = null;
	}
	
	public synchronized void makeFixedAtEnd(double end, double variance) {
		type = FIXED_AT_END;
		this.end = end;
		this.variance = variance;
		endDimensions = 1;
		pathFunction = null;
	}
	
	public synchronized void makeFixed(double start, double end, double variance) {
		type = FIXED;
		this.start = start;
		this.end = end;
		this.variance = variance;
		endDimensions = 0;
		pathFunction = null;
	}
	
	public synchronized void makeForTrace(double variance) {
		type = TRACE;
		this.variance = variance;
		endDimensions = 1;
		pathFunction = null;
	}
	
	@Override
	public synchronized int inputDimension() {
		return steps - 1 + endDimensions;
	}
	
	@Override
	public synchronized int outputDimension() {
		return steps + 1;
	}
	
	@Override
	public synchronized double[] evaluate(double[] x) {
		checkArgument(x);
		
		double[] path;
		
		if (pathFunction == null) {
			makePathFunction();
		}
		
		if (type == TRACE) {
			double[] xx = new double[x.length - 1];
			System.arraycopy(x, 1, xx, 0, xx.length);
			path = pathFunction.evaluate(xx);
			shiftPath(path, x[0], x[0]);
		} else {
			path = pathFunction.evaluate(x);
			if (type == FIXED_AT_START) {
				shiftPath(path, start, start);
			} else if (type == FIXED_AT_END) {
				path = Util.reverseArray(path);
				shiftPath(path, end, end);
			} else if (type == FIXED) {
				double[] xx = new double[x.length];
				System.arraycopy(x, 0, xx, 0, xx.length);
				shiftPath(path, start, end);
			} else {
				throw new IllegalStateException("Illegal type");
			}
		}
		
		return path;
	}
	
	private void shiftPath(double[] path, double start, double end) {
		double slope = (end - start) / steps;
		double shift = start;
		
		for (int i = 0; i < path.length; i++) {
			path[i] += shift;
			shift += slope;
		}
	}
	
	private void makePathFunction() {
		boolean open = type == FIXED_AT_START || type == FIXED_AT_END;
		
		if (generatorType == BRIDGE) {
			pathFunction = new BrownianBridge(variance, steps, open);
		} else if (generatorType == WALK) {
			if (open) {
				pathFunction = new BrownianWalk(variance, steps);
			} else {
				throw new IllegalArgumentException("Type WALK is not compatible with a closed path");
			}
		} else {
			throw new IllegalArgumentException("Invalid type");
		}
	}
	
	@Override
	public synchronized String toString() {
		return "Brownian Path(type=" + type + ",steps=" + steps + ",v=" + variance + ",x=" + start + ",y=" + end + ")";
	}
}
