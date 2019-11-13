package de.labathome;

import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

import static de.labathome.Profiles.AtanProfile;
import java.util.Locale;

public class TestProfiles {

	@Test
	public final void testAtanProfileAnalyticIntegralAndDerivative() {

		double am0 = 13.0;
		
		double minA = 0.5;
		double maxA = 5.0;
		int numA = 11;
		
		double minB = 0.5;
		double maxB = 5.0;
		int numB = 9;
		
		int N = 20;
		
		double cubatureRelTol = 1.0e-12;
		double cubatureAbsTol = Double.NaN;
		double checkTol = 1.0e-9;

		double[] aRange = new double[numA];
		for (int i = 0; i < numA; ++i) {
			aRange[i] = minA + i * (maxA - minA) / (numA - 1.0);
		}

		double[] bRange = new double[numB];
		for (int i = 0; i < numB; ++i) {
			bRange[i] = minB + i * (maxB - minB) / (numB - 1.0);
		}

		double[] testRange = new double[N];
		for (int i = 0; i < N; ++i) {
			testRange[i] = 0.01 + i * (1.0 - 0.01) / (N - 1.0);
		}

		// use Cubature for adaptive integration as reference; needs multi-dimensional
		// wrapper around one-dimensional AtanProfile integrand
		class AtanIntegrand implements Cubature.Integrand {
			@Override
			public double[][] eval(double[][] x, Object fdata) {
				if (fdata instanceof AtanProfile) {
					return new double[][] { ((AtanProfile) fdata).eval(x[0]) };
				} else {
					throw new RuntimeException("fdata has to be of type AtanProfile");
				}
			}
		}
		AtanIntegrand atanIntegrand = new AtanIntegrand();

		for (int iA = 0; iA < numA; ++iA) {
			for (int iB = 0; iB < numB; ++iB) {
			
				// System.out.println("a="+aRange[iA]+" b="+bRange[iB]);

				AtanProfile p = new AtanProfile(am0, aRange[iA], bRange[iB]);

				double[] numerical = new double[N];
				double[] analytical = new double[N];

				// loop over lower boundary; does not range up to 1
				for (int i = 0; i < N-1; ++i) {
					double s0 = testRange[i];
					
					// loop over upper boundary; has to be greater than or equal to lower boundary
					for (int j = i+1; j < N; ++j) {
						double s1 = testRange[j];
	
						analytical[i] = p.integral(s0, s1);
	
						double[][] valErr = Cubature.integrate(atanIntegrand, "eval", new double[] { s0 }, new double[] { s1 },
								cubatureRelTol, cubatureAbsTol, Cubature.Error.INDIVIDUAL, 1000000, p);
						numerical[i] = valErr[0][0];
						
						double relDev = (analytical[i] - numerical[i])/numerical[i];
						
						boolean success = Math.abs(relDev) < checkTol;
						if (!success) {
							System.out.println("failed to satisfy tolerance requirement");
							System.out.println(" am0 = "+am0);
							System.out.println("  iA = "+iA);
							System.out.println("   a = "+aRange[iA]);
							System.out.println("  iB = "+iB);
							System.out.println("   b = "+bRange[iB]);
							System.out.println("   i = "+i);
							System.out.println("  s0 = "+s0);
							System.out.println("   j = "+j);
							System.out.println("  s1 = "+s1);
							System.out.println("  rel. dev. = "+relDev);
						}
						assertTrue(success);
						
						System.out.println(String.format(Locale.ENGLISH,
								"integral from %.3e to %3e : analytical = %.3e ; numerical = %.3e ; rel. dev = %.3e", s0,
								s1, analytical[i], numerical[i], relDev));
					}
				}
			}
		}
	}
}
