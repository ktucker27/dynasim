package exe;

import handlers.RingSteadyStateTest;
import handlers.SynchSteadyStateTerminator;
import handlers.WriteHandler;

import java.io.FileNotFoundException;
import java.io.UnsupportedEncodingException;
import java.util.Random;

import ode.SynchODEs;

import org.apache.commons.math3.ode.nonstiff.AdamsMoultonIntegrator;

import coupling.RingCoupling;

public class SimSynchSingle {

    /**
     * @param args
     * @throws UnsupportedEncodingException 
     * @throws FileNotFoundException 
     */
    public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException {
        long startTime = System.nanoTime();

        int n = 103;
        double h = 0.0005;
        double a = 0.025;
        
        double[] d = new double[n];
        for(int i = 0; i < n; ++i) {
            d[i] = 0.0;
        }
        
        double[] y0 = new double[3*n];
        double[] y = new double[3*n];
        
        Random rand = new Random(5);
        for(int i = 0; i < n; ++i) {
            y0[i] = rand.nextDouble();
            y0[n+i] = rand.nextDouble();
            y0[2*n+i] = rand.nextDouble()*2*Math.PI - Math.PI;
        }

        double t0 = 0.0;
        double t = 10.0;
        
        double w = 10.0;

        WriteHandler writeHandler = new WriteHandler("/Users/kristophertucker/Google Drive/Research/Synch/output/r_vs_t.txt", new int[] {n+n/2});
        
        RingSteadyStateTest test = new RingSteadyStateTest();
        test.setQuietMode(false);
        SynchSteadyStateTerminator term = new SynchSteadyStateTerminator(test, 0.5, 500, 1500000);
        
        AdamsMoultonIntegrator integrator = new AdamsMoultonIntegrator(2, 1.0e-18, h, 1.0e-3, 1.0e-2);
        integrator.addStepHandler(writeHandler);
        integrator.addStepHandler(term.getDetector());
        integrator.addEventHandler(term, Double.POSITIVE_INFINITY, 1.0e-12, 100);
        
        RingCoupling coupling = new RingCoupling(n,a);
        SynchODEs odes = new SynchODEs(n, 1, w, coupling, d);
        try {
            integrator.integrate(odes, t0, y0, t, y);
            System.out.println("R:");
            for(int i = n; i < 2*n; ++i) {
                System.out.print(y[i] + " ");
            }
            System.out.print("\n");
        } catch(Exception ex) {
            System.out.println(ex.getMessage());
        }

        long endTime = System.nanoTime();

        System.out.println("Run time: " + (endTime - startTime)/1.0e9 + " seconds");
    }
}
