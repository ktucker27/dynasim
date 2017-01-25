package exe;

import integrator.IntegratorFactory;
import integrator.ThreadPoolIntegrator;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.Random;

import ode.SynchMeanFieldODEs;

import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.AdamsBashforthIntegrator;

import utils.OrderParameterSolution;

public class SimSynchMeanField {

    /**
     * @param args
     * @throws UnsupportedEncodingException 
     * @throws FileNotFoundException 
     */
    public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException {
        long startTime = System.nanoTime();

        int n = 70;
        //double h = 0.01;
        
        double[] d = new double[n];
        for(int i = 0; i < n; ++i) {
            d[i] = 0.0;
        }
        
        double[] y0 = new double[3*n];
        
        Random rand = new Random(1);
        for(int i = 0; i < n; ++i) {
            y0[i] = rand.nextDouble();
            y0[n+i] = rand.nextDouble();
            y0[2*n+i] = rand.nextDouble()*2*Math.PI - Math.PI;
        }

        double t0 = 0.0;
        double t = 6.0;
        
        double hy = 69.0/200.0;
        double hx = 39.0/200.0;
        
        IntegratorFactory factory = new IntegratorFactory() {
            double myMinStep = 1.0e-18;
            double myMaxStep = 0.01;
            
            @Override
            public FirstOrderIntegrator newIntegrator() {
                return new AdamsBashforthIntegrator(2, myMinStep, myMaxStep, 1.0e-3, 1.0e-2);
            }
        };

        ThreadPoolIntegrator integrator = new ThreadPoolIntegrator(8, factory);
        integrator.setQuietMode(false);

        OrderParameterSolution[][] r = new OrderParameterSolution[201][201];
        for(int i = 0; i < 201; ++i) {
            for(int j = 0; j < 201; ++j) {
                r[i][j] = new OrderParameterSolution();
                SynchMeanFieldODEs odes = new SynchMeanFieldODEs(n, 1, 1.0 + j*hx, 1.0 + i*hy, 0.0, d);
                integrator.addIvp(odes, t0, y0, t, r[i][j]);
            }
        }

        boolean success = integrator.waitForFinished(3600);

        long endTime = System.nanoTime();

        if(success) {
            PrintWriter writer = new PrintWriter("/Users/kristophertucker/Google Drive/Research/Synch/output/mesh.txt", "UTF-8");
            for(int i = 0; i < 201; ++i) {
                for(int j = 0; j < 201; ++j) {
                    writer.print(r[i][j].getOrderParam() + " ");
                }
                writer.print("\n");
            }
            writer.close();
        }

        System.out.println("Run time: " + (endTime - startTime)/1.0e9 + " seconds");
    }
}
