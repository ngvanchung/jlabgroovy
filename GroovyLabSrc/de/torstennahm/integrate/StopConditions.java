/*
 * Created on May 3, 2004
 */
package de.torstennahm.integrate;

import java.util.Iterator;
import java.util.List;

/**
 * This class provides several common stopping conditions as
 * static subclasses.
 * 
 * This class cannot be instantiated.
 * 
 * @author Torsten Nahm
 */
public class StopConditions {
	private StopConditions() {}
	
	/**
	 * Signals the integration should be stopped if the specified timeout is reached.
	 * 
	 * @author Torsten Nahm
	 */
	static public class UntilTimeOut implements StopCondition {
		private long timeOut;
	
		/**
		 * Constructs the stop condition.
		 * Time is counted from the call of this constructor.
		 * 
		 * @param ms timeout for integration
		 */
		public UntilTimeOut(long ms) {
			timeOut = System.currentTimeMillis() + ms;
		}
	
		public boolean stop(IntegrationResult result) {
			return System.currentTimeMillis() >= timeOut;
		}
	
		public String getConditionString() {
			return "*** Timeout ***";
		}
	}

	/**
	 * Signals the integration should be stopped if the integrand functions
	 * has been called at least the specified number of times.
	 */ 
	static public class UntilCallsReached implements StopCondition {
		private long minCalls; 

		/**
		 * Construct the stop condition.
		 * 
		 * @param calls minimum number of evaluations
		 */
		public UntilCallsReached(long calls) {
			this.minCalls = calls;
		}
	
		public boolean stop(IntegrationResult result) {
			return result.functionCalls() >= minCalls;
		}
	
		public String getConditionString() {
			return "Function evaluations reached";
		}
	}
	
	/**
	 * Signals the integration should be stopped if the specified absolute
	 * tolerance has been reached.
	 * 
	 * That is, this class checks whether the error estimate has an absolute value
	 * equal to or less than the given tolerance.
	 * <p>
	 * If the check fails, the current number of function evaluations <i>n</i>
	 * will be stored. The next check will only be performed when the
	 * number of function evaluations is at least <i>n*1.5</i>. This is
	 * an effort to mitigate stopping bias: since an integration that is
	 * stopped will not be restarted, but an integration that is not
	 * stopped will be checked for stopping again, and since the
	 * error is an estimate underlying random fluctuations, the
	 * net effect is a premature stopping, which becomes worse
	 * the more often the error is sampled.
	 */
	static public class UntilAbsTol implements StopCondition {
		private double absTol;
		private long nextCheck = 10;
		
		/**
		 * Construct the stop condition.
		 * 
		 * @param absoluteTolerance tolerance to be reached
		 */
		public UntilAbsTol(double absoluteTolerance) {
			absTol = absoluteTolerance;
		}
	
		public boolean stop(IntegrationResult result) {
			if (result.functionCalls() >= nextCheck) {
				if (result.errorEstimate() <= absTol) {
					return true;
				}
				nextCheck = (result.functionCalls() * 3 / 2) + 1;
			}
			
			return false;
		}
	
		public String getConditionString() {
			return "Absolute Tolerance achieved";
		}
	}
	
	/**
	 * Signals the integration should be stopped if the specified absolute
	 * tolerance has been reached.
	 * 
	 * That is, this class checks whether the error estimate divided by the estimate
	 * for the value of the integral has an absolute value less than
	 * the given tolerance.
	 * <p>
	 * If the check fails, the current number of function evaluations <i>n</i>
	 * will be stored. The next check will only be performed when the
	 * number of function evaluations is at least <i>n*1.5</i>. This is
	 * an effort to mitigate stopping bias: since an integration that is
	 * stopped will not be restarted, but an integration that is not
	 * stopped will be checked for stopping again, and since the
	 * error is an estimate underlying random fluctuations, the
	 * net effect is a premature stopping, which becomes worse
	 * the more often the error is sampled.
	 */
	static public class UntilRelTol implements StopCondition {
		private double relTol;

		/**
		 * Construct the stop condition.
		 * 
		 * @param relativeTolerance tolerance to be reached
		 */
		public UntilRelTol(double relativeTolerance) {
			relTol = relativeTolerance;
		}

		public boolean stop(IntegrationResult result) {
			return (result.errorEstimate() / Math.abs(result.value())) <= relTol;
		}
	
		public String getConditionString() {
			return "Relative Tolerance achieved";
		}
	}
	
	/**
	 * Checks whether any of a list of conditions has been fulfilled.
	 */
	static public class MultiStopCondition implements StopCondition {
		private List<StopCondition> conditions;
		private StopCondition stopper = null;
		
		/**
		 * Constructs the stop condition.
		 * 
		 * @param conditions list of conditions
		 */
		public MultiStopCondition(List<StopCondition> conditions) {
			this.conditions = conditions;
		}

		public boolean stop(IntegrationResult result) {
			for (Iterator<StopCondition> iter = conditions.iterator(); iter.hasNext();) {
				StopCondition ci = iter.next();
				if (ci.stop(result)) {
					stopper = ci;
					return true;
				}
			}
			
			return false;
		}
		
		/**
		 * Returns which condition has caused the <code>stop</code> method to return
		 * <code>true</code>. If <code>stop</code> has not yet returned <code>true</code>,
		 * <code>null</code> will be returned.
		 * 
		 * @return condition that caused the stop, or <code>null</code>
		 */
		public StopCondition getStopper() {
			return stopper;
		}
	
		public String getConditionString() {
			StringBuffer sb = new StringBuffer();
			for (Iterator<StopCondition> iter = conditions.iterator(); iter.hasNext();) {
				sb.append((iter.next()).getConditionString());
				if (iter.hasNext()) {
					sb.append(" or ");
				}
			}
			
			return sb.toString();
		}
	}
}