package exe;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.Random;

import ode.SynchMeanFieldODEs;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.ode.nonstiff.EulerIntegrator;

public class SimSynchMeanFieldST {

    /**
     * @param args
     * @throws UnsupportedEncodingException 
     * @throws FileNotFoundException 
     */
    public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException {
        long startTime = System.nanoTime();

        double hy = 69.0/200.0;
        double hx = 39.0/200.0;

        double[][] r = new double[201][201];
        for(int i = 0; i < 201; ++i) {
            for(int j = 0; j < 201; ++j) {
                r[i][j] = runSim(1.0 + j*hx, 1.0 + i*hy, 0.0);
            }
            System.out.println("Finished row " + i);
        }

        long endTime = System.nanoTime();

        PrintWriter writer = new PrintWriter("/Users/kristophertucker/output/mesh.txt", "UTF-8");
        for(int i = 0; i < 201; ++i) {
            for(int j = 0; j < 201; ++j) {
                writer.print(r[i][j] + " ");
            }
            writer.print("\n");
        }
        writer.close();

        System.out.println("Run time: " + (endTime - startTime)/1.0e9 + " seconds");
    }

    public static double runSim(double w, double feff, double geff) {
        int n = 70;
        double h = 0.01;
        double[] d = new double[n];
        for(int i = 0; i < n; ++i) {
            d[i] = 0.0;
        }
        SynchMeanFieldODEs odes = new SynchMeanFieldODEs(n, 1, w, feff, geff, d);
        
        EulerIntegrator integrator = new EulerIntegrator(h);
        
        double[] y = new double[3*n];
        double[] y0 = new double[3*n];
        
        Random rand = new Random();
        for(int i = 0; i < n; ++i) {
            y0[i] = rand.nextDouble();
            y0[n+i] = rand.nextDouble();
            y0[2*n+i] = rand.nextDouble()*2*Math.PI - Math.PI;
        }
        
        double t = 6.0;
        integrator.integrate(odes, 0, y0, t, y);
        
        // Compute the order parameter
        Complex z = new Complex(0.0);
        for(int i = 0; i < n; ++i) {
            z = z.add(Complex.I.multiply(y[2*n+i]).exp().multiply(y[n+i]));
        }
        z = z.multiply(1/(double)n);

        return z.abs();
        //double psi = z.getArgument();
    }
}
