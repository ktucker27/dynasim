package exe;

import integrator.EulerIntegratorFactory;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.Random;

import ode.SynchMeanFieldODEs;

import org.apache.commons.math3.complex.Complex;

import utils.ThreadPoolIntegrator;

public class SimSynchMeanField {

    /**
     * @param args
     * @throws UnsupportedEncodingException 
     * @throws FileNotFoundException 
     */
    public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException {
        long startTime = System.nanoTime();

        int n = 70;
        double h = 0.01;
        
        double[] d = new double[n];
        for(int i = 0; i < n; ++i) {
            d[i] = 0.0;
        }
        
        double[] y0 = new double[3*n];
        
        Random rand = new Random();
        for(int i = 0; i < n; ++i) {
            y0[i] = rand.nextDouble();
            y0[n+i] = rand.nextDouble();
            y0[2*n+i] = rand.nextDouble()*2*Math.PI - Math.PI;
        }

        double t0 = 0.0;
        double t = 6.0;
        
        double hy = 69.0/200.0;
        double hx = 39.0/200.0;

        ThreadPoolIntegrator integrator = new ThreadPoolIntegrator(8, new EulerIntegratorFactory(h));
        integrator.setQuietMode(false);

        double[][] r = new double[201][201];
        double[][][] y = new double[201][201][3*n];
        for(int i = 0; i < 201; ++i) {
            for(int j = 0; j < 201; ++j) {
                SynchMeanFieldODEs odes = new SynchMeanFieldODEs(n, 1, 1.0 + j*hx, 1.0 + i*hy, d);
                integrator.addIvp(odes, t0, y0, t, y[i][j]);
            }
        }

        boolean success = integrator.waitForFinished(3600);

        if(success) {
            // Compute order parameters
            for(int i = 0; i < 201; ++i) {
                for(int j = 0; j < 201; ++j) {
                    r[i][j] = compOrderParam(y[i][j]);
                }
            }
        }
        
        long endTime = System.nanoTime();

        if(success) {
            PrintWriter writer = new PrintWriter("/Users/kristophertucker/output/mesh_thread.txt", "UTF-8");
            for(int i = 0; i < 201; ++i) {
                for(int j = 0; j < 201; ++j) {
                    writer.print(r[i][j] + " ");
                }
                writer.print("\n");
            }
            writer.close();
        }

        System.out.println("Run time: " + (endTime - startTime)/1.0e9 + " seconds");
    }

    public static double compOrderParam(double[] y) {
        // Compute the order parameter
        int n = y.length/3;
        Complex z = new Complex(0.0);
        for(int i = 0; i < n; ++i) {
            z = z.add(Complex.I.multiply(y[2*n+i]).exp().multiply(y[n+i]));
        }
        z = z.multiply(1/(double)n);

        return z.abs();
        //double psi = z.getArgument();
    }
}
