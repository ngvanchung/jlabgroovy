/*
 * Created on Feb 3, 2004
 */
package de.torstennahm.integrate;

/**
 * This inteface provides a mechanism for deciding when integration should stop.
 *
 * @author Torsten Nahm
 */
	
public interface StopCondition {
	/**
	 * Returns true if the integration should be stopped because of this condition.
	 * 
	 * @param result integration result
	 * @return true if condition is fulfilled
	 */
	boolean stop(IntegrationResult result);
	
	/**
	 * Returns a string describing the fulfillment of the condition. For example,
	 * if stopping when the error tolerance is below the requirement, the string
	 * might be "Error tolerance reached".
	 * 
	 * @return string describing condition when it is fulfilled
	 */
	String getConditionString();
}