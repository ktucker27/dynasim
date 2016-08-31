package exe;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.Random;

import ode.SynchODEs;

import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.nonstiff.AdamsMoultonIntegrator;
import org.apache.commons.math3.ode.sampling.StepHandler;
import org.apache.commons.math3.ode.sampling.StepInterpolator;

import utils.SynchSteadyStateTerminator;

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
        
        Random rand = new Random(2);
        for(int i = 0; i < n; ++i) {
            y0[i] = rand.nextDouble();
            y0[n+i] = rand.nextDouble();
            y0[2*n+i] = rand.nextDouble()*2*Math.PI - Math.PI;
        }

        double t0 = 0.0;
        double t = 10.0;
        
        double w = 10.0;

        StepHandler writeHandler = new StepHandler() {
            PrintWriter writer = new PrintWriter("/Users/kristophertucker/Google Drive/Research/Synch/output/r_vs_t.txt", "UTF-8");
            
            @Override
            public void init(double t0, double[] y0, double t) {
                int n = y0.length/3;
                writer.print("0, " + y0[n] + "\n");
            }
            
            @Override
            public void handleStep(StepInterpolator interpolator, boolean isLast)
                    throws MaxCountExceededException {
                double[] y = interpolator.getInterpolatedState();
                int n = y.length/3;
                
                writer.print(interpolator.getInterpolatedTime() + ", ");
                writer.print(y[n]);
                writer.print("\n");
                
                if(!isLast) {
                    writer.flush();
                } else {
                    writer.close();
                }
            }
        };
        
        SynchSteadyStateTerminator term = new SynchSteadyStateTerminator(500);
        
        AdamsMoultonIntegrator integrator = new AdamsMoultonIntegrator(2, 1.0e-18, h, 1.0e-3, 1.0e-2);
        integrator.addStepHandler(writeHandler);
        integrator.addStepHandler(term.getDetector());
        integrator.addEventHandler(term, Double.POSITIVE_INFINITY, 1.0e-12, 100);
        
        RingCoupling coupling = new RingCoupling(n,a);
        SynchODEs odes = new SynchODEs(n, 1, w, coupling, d);
        try {
            integrator.integrate(odes, t0, y0, t, y);
        } catch(Exception ex) {
            System.out.println(ex.getMessage());
        }

        long endTime = System.nanoTime();

        System.out.println("Run time: " + (endTime - startTime)/1.0e9 + " seconds");
    }
}
