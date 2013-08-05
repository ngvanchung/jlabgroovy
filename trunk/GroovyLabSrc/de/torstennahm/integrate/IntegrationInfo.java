/*
 * Created on Aug 4, 2004
 */
package de.torstennahm.integrate;

/**
 * Provides additional information for an <code>IntegrationResult</code>.
 * 
 * This class provides only a human-readable text description of the information.
 * Other information should be provided by subclasses with additional methods.
 * <p>
 * Implementor's note: The integration info objects are stored in a set. In this
 * way, if the same information is generated several times, and added to the set,
 * only one copy is retained. This prevents the used from getting the same
 * message potentially hundreds of times. For this to work, only the same information
 * must indeed compare as equal in the set. Subclasses therefore need to
 * override <code>equals</code>.
 * 
 * @author Torsten Nahm
 */
public class IntegrationInfo {
	protected String errorDescription;
	
	public IntegrationInfo(String textInfo) {
		this.errorDescription = textInfo;
	}
	
	public String textInfo() {
		return errorDescription;
	}
	
	@Override
	public boolean equals(Object o) {
		return (o.getClass() == IntegrationInfo.class) && textInfo().equals(((IntegrationInfo) o).textInfo());
	}
	
	@Override
	public int hashCode() {
		return errorDescription.hashCode();
	}
	
	@Override
	public String toString() {
		return textInfo();
	}
}
