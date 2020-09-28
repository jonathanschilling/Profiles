package de.labathome;

/**
 * This class presents bundles functions occuring for e.g. radial profiles of plasma parameters in a nuclear fusion plasma
 * or any other case where some quantitity is representable by a continous, analytical profile over a domain [0,1].
 * 
 * @author Jonathan Schilling (jonathan.schilling@mail.de)
 */
public class Profiles {
	
	/** A parameterized profile which two shaping parameters. */
	public static class AtanProfile {
		
		private double am0, a, b, p2, p3;
		
		/**
		 * Initialize a new instance of the parameterized profile.
		 * @param am0 central value, i.e. the one returned for for x=0
		 * @param a steepness parameter
		 * @param b location of largest gradient
		 */
		public AtanProfile(double am0, double a, double b) {
			this.am0 = am0; this.a = a; this.b = b;
			this.p2 = Math.atan(a*b-a);
			this.p3 = Math.atan(a*b);
		}
		
		/**
		 * Evaluate the profile function at a number of locations.
		 * Tested for x \in [0,1]
		 * @param x [N] locations at which to evaluate the profile function
		 * @return f(x) [N] the profile evaluated at locations {@code x}
		 */
		public double[] eval(double[] x) {
			double[] p = new double[x.length];
			double p1;
			for (int i=0; i<x.length; ++i) {
				p1 = Math.atan(a*(b-x[i]));
				p[i] = am0*(p1 - p2)/(p3 - p2);
			}
			return p;
		}
		
		/**
		 * Analytical integral of the profile function.
		 * Tested for x0 &lt; x1 where (x0,x1) \in [0,1]^2
		 * @param x0 lower integration boundary
		 * @param x1 upper integration boundary
		 * @return \int_x0^x1 f(x) dx
		 */
		public double integral(double x0, double x1) {
			double ps  = (x0-x1)*p2;
			double ps0 = (b-x0)*Math.atan(a*(b-x0));
			double ps1 = (b-x1)*Math.atan(a*(b-x1));
			double ls0 = Math.log(1.0+a*a*(b-x0)*(b-x0))/(2.0*a);
			double ls1 = Math.log(1.0+a*a*(b-x1)*(b-x1))/(2.0*a);
			return am0*(ps + ps0 - ps1 - ls0 + ls1)/(p3 - p2);
		}
		
		/**
		 * Analytical derivative of the profile function.
		 * Untested.
		 * @param x location at which to compute the derivative.
		 * @return df/dx at x
		 */
		public double derivative(double x) {
			double denom = (1.0 + a*a*(b-x)*(b-x))*(p3 - p2);
			return -a*am0/denom;
		}
	};
}
