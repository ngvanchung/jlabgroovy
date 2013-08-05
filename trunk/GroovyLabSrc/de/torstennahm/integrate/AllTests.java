/*
 * Created on Nov 6, 2003
 */
package de.torstennahm.integrate;

import junit.framework.Test;
import junit.framework.TestSuite;

/**
 * JUnit suite to run all JUnit tests in de.torstennahm.integrate and subpackages.
 * 
 * @author Torsten Nahm
 */
public class AllTests {
	/**
	 * Create test suite.
	 * 
	 * @return test suite
	 */
	public static Test suite() {
		TestSuite suite =
			new TestSuite("Test for de.torstennahm.integrate");
		//$JUnit-BEGIN$
		suite.addTest(new TestSuite(de.torstennahm.integrate.quadratureformula.JUnitTest.class));
		suite.addTest(new TestSuite(de.torstennahm.integrate.sparse.index.JUnitIndexTest.class));
		//$JUnit-END$
		return suite;
	}
}
