/*
 * Created on Jul 17, 2003
 */
package de.torstennahm.integrate;

/**
 * This exception is thrown if the integration cannot be completed.
 * 
 * @author Torsten Nahm
 */
public class IntegrationFailedException extends Exception {
    static final long serialVersionUID = 1094686409050309601L;
    
	/**
	 * Constructs a new exception with <code>null</code> as its detail message.
	 */
	public IntegrationFailedException() {
		super();
	}
	
	/**
	 * Constructs a new exception with the specified detail message.
	 * 
	 * @param string detail message
	 */
	public IntegrationFailedException(String string) {
		super(string);
	}
}
