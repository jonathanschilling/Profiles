package examples;

import java.util.Locale;

import de.labathome.Profiles.AtanProfile;

public class Examples {

	public static void main(String[] args) {
		ex_AtanProfile();
	}
	
	public static void ex_AtanProfile() {
		
		double am0 = 120;
		double a = 10;
		double b = 0.2;
		
		AtanProfile p = new AtanProfile(am0, a, b);
		
		int N = 20;
		double[] x = new double[N];
		for (int i=0; i<N; ++i) {
			x[i] = i/(N-1.0);
		}
		
		double[] y = p.eval(x);
		
		for (int i=0; i<N; ++i) {
			System.out.println(String.format(Locale.ENGLISH,
					"%.6e %.6e", x[i], y[i]
			));
		}
		
		// put the output into a text file, e.g. "p.dat" and plot it using e.g. gnuplot:
		// gnuplot> plot 'p.dat' u 1:2 w lp
	}

}
