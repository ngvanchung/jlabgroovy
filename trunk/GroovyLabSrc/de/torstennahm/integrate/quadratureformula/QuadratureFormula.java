/*
 * Created on Jul 6, 2003
 */
package de.torstennahm.integrate.quadratureformula;

import java.util.Iterator;
import java.util.List;

/**
 * This class represents a quadrature formula, given by a list of nodes and weights.
 * 
 * The nodes are sorted in ascending order, so lower indices correspond to
 * numerically smaller nodes.
 * <p>
 * According to general <code>Generator</code> contract, this class is thread-safe.
 * 
 * @author Torsten Nahm
 */
public class QuadratureFormula {
	final private double[] nodes;
	final private double[] weights;

	/**
	 * Constructs the quadrature formula from arrays. The array of nodes must be sorted in
	 * ascending order.
	 * 
	 * @param nodes array of nodes
	 * @param weights array of weights
	 * @throws IllegalArgumentException if the arrays have different lengths
	 * @throws IllegalArgumentException if the nodes are not sorted
	 */
	public QuadratureFormula(double[] nodes, double[] weights) {
		this.nodes = nodes.clone();
		this.weights = weights.clone();
		checkIntegrity();
	}
	
	/**
	 * Constructs the quadrature formula from lists of <code>Double</code>s. The list of nodes must be
	 * sorted in ascending order.
	 * 
	 * @param nodeList list of nodes
	 * @param weightList list of weights
	 * @throws IllegalArgumentException if the lists have different lengths
	 * @throws IllegalArgumentException if the nodes are not sorted
	 */
	public QuadratureFormula(List<Double> nodeList, List<Double> weightList) {
		nodes = new double[nodeList.size()];
		weights = new double[weightList.size()];
		int i = 0;
		for (Iterator<Double> iter = nodeList.iterator(); iter.hasNext(); ) {
			nodes[i++] = iter.next().doubleValue();
		}
		i = 0;
		for (Iterator<Double> iter = weightList.iterator(); iter.hasNext(); ) {
			weights[i++] = iter.next().doubleValue();
		}
		checkIntegrity();
	}
	
	/**
	 * Check if the nodes and weights the object is constructed from
	 * have the same length, and if the nodes are stored.
	 */
	private void checkIntegrity() {
		if (nodes.length != weights.length) {
			throw new IllegalArgumentException("Node and weight must have the same length");
		}
	
		double curNode = Double.NEGATIVE_INFINITY;
		for (int i = 0; i < nodes.length; i++) {
			if (nodes[i] < curNode) {
				throw new IllegalArgumentException("Nodes must be sorted");
			}
			curNode = nodes[i];
		}
	}

	/**
	 * Returns number of (node, weight) pairs of the formula.
	 * 
	 * @return number of nodes in the forumla
	 */
	public int getSize() {
		return nodes.length;
	}

	/**
	 * Returns the node of the specified (node, weight) pair.
	 * 
	 * @param index of the pair
	 * @return node with the given index
	 */
	public double getNode(int index) {
		return nodes[index];
	}

	/**
	 * Returns the weight of the specified (node, weight) pair.
	 * 
	 * @param index of the pair
	 * @return weight with the given index
	 */
	public double getWeight(int index) {
		return weights[index];
	}
	
	/**
	 * Returns the nodes as an array.
	 * 
	 * @return nodes array
	 */
	public double[] getNodesArray() {
		return nodes.clone();
	}
	
	/**
	 * Returns the weights as an array.
	 * 
	 * @return weights array
	 */
	public double[] getWeightsArray() {
		return weights.clone();
	}
}
