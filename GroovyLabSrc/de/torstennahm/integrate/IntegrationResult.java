/*
 * Created on Jun 24, 2004
 */
package de.torstennahm.integrate;

import java.util.Set;

/**
 * This interface models the result of an integration.
 * 
 * The result consists of the estimated value,
 * its estimated error and the number of function evaluations used for the integration.
 * 
 * @author Torsten Nahm
 */
public interface IntegrationResult {
	/**
	 * Returns the estimated integral value. <code>Double.NaN</code> is returned if
	 * no estimate is available.
	 *  
	 * @return integral value
	 */
	double value();
	
	/**
	 * Returns the error estimate for the integral value. <code>Double.NaN</code> is returned
	 * if no estimate is available.
	 * 
	 * @return error estimate of integral value
	 */
	double errorEstimate();
	
	/**
	 * Returns the number of times the integrand function has been evaluated.
	 * 
	 * @return number of evaluations
	 */
	long functionCalls();
	
	/**
	 * Returns a set of integration information that
	 * provide supplementary or additional information on the integration
	 * process and result.
	 * 
	 * The set may be empty.
	 * 
	 * @return set of integration information
	 */
	Set<IntegrationInfo> supplementalInfo();
}
