/*
 * Created on Aug 5, 2004
 */
package de.torstennahm.integrate.sparse;

import de.torstennahm.integrate.IntegrationInfo;
import de.torstennahm.integrate.sparse.evaluateindex.Evaluator;
import de.torstennahm.integrate.sparse.index.Index;

/**
 * Signifies that an integration index could not be evaluated during the integration process.
 * 
 * @author Torsten Nahm
 */
public class CouldNotEvaluateInfo extends IntegrationInfo {
	public CouldNotEvaluateInfo(Evaluator evaluator, Index index) {
		super("Could not evaluate an integration index. This error first occurred for the index " + index + ".\n" +
			  "Although integration could be completed, this condition implies slower convergence " +
			  "and error overestimation.\n" +
			  "Please check the index evaluator for possible causes.");
	}
	
	@Override
	public boolean equals(Object o) {
		return o.getClass() == CouldNotEvaluateInfo.class;
	}
	
	@Override
	public int hashCode() {
		return CouldNotEvaluateInfo.class.hashCode();
	}
}
