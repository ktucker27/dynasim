package exe;

import handlers.MeanFieldSteadyStateTest;
import handlers.SynchSteadyStateTerminator;
import handlers.WriteHandlerOP;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.Random;

import ode.SynchMeanFieldODEs;

import org.apache.commons.math3.distribution.TDistribution;
import org.apache.commons.math3.ode.nonstiff.AdamsMoultonIntegrator;
import org.apache.commons.math3.random.Well19937c;

import utils.OrderParameterSolution;

public class SimBistability {

    /**
     * @param args
     * @throws UnsupportedEncodingException 
     * @throws FileNotFoundException 
     */
    public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException {
        long startTime = System.nanoTime();

        // Simulation parameters
        int n = 103;
        double h = 0.001;
        double delta = 2.0;
        double w = 2.0;
        double gamma = 1.0;
        
        double fmin = 17.4;
        double fmax = 15.0;
        int numf = 35;
        
        double deltaf = (fmax - fmin)/(double)(numf - 1);
        
        double geff = 2.0*Math.sqrt((gamma+w+delta)*Math.pow(gamma+w,2)*Math.pow(gamma+2*delta+w, 2)/(delta*gamma*gamma*Math.pow(w-gamma, 2)));
        System.out.println("g_eff: " + geff);
        
        // Get natural frequencies from Cauchy distribution
        TDistribution tdist = new TDistribution(new Well19937c(1), 1);
        double[] d = tdist.sample(n);
        for(int i = 0; i < n; ++i) {
            d[i] = d[i]*delta;
        }
        
        double[] y0 = new double[3*n];
        double[] y = new double[3*n];
        
        Random rand = new Random(1);
        for(int i = 0; i < n; ++i) {
            y0[i] = rand.nextDouble();
            y0[n+i] = rand.nextDouble();
            y0[2*n+i] = rand.nextDouble()*2*Math.PI - Math.PI;
        }

        double t0 = 0.0;
        double t = 1000.0;

        PrintWriter writer = new PrintWriter("/Users/kristophertucker/Google Drive/Research/Synch/output/bistability.txt", "UTF-8");
        for(int i = 0; i < numf; ++i) {
            MeanFieldSteadyStateTest test = new MeanFieldSteadyStateTest(1.0e-5);
            SynchSteadyStateTerminator term = new SynchSteadyStateTerminator(test, 3.0, 200, 5000000);
            AdamsMoultonIntegrator integrator = new AdamsMoultonIntegrator(2, 1.0e-18, h, 1.0e-3, 1.0e-2);
            integrator.addStepHandler(term.getDetector());
            integrator.addEventHandler(term, Double.POSITIVE_INFINITY, 1.0e-12, 100);
            
            double feff = fmin + i*deltaf;
            SynchMeanFieldODEs odes = new SynchMeanFieldODEs(n, gamma, w, feff, geff, d);

            if(feff > 16.46 && feff < 16.5) {
                System.out.println("Recording for feff = " + feff);
                WriteHandlerOP writeHandler = new WriteHandlerOP("/Users/kristophertucker/Google Drive/Research/Synch/output/z_vs_t.txt");
                integrator.addStepHandler(writeHandler);
            }
            
            try {
                System.out.println(i + " " + feff);
                integrator.integrate(odes, t0, y0, t, y);
                System.out.println(OrderParameterSolution.compOrderParam(y));
                
                writer.println(feff + ", " + OrderParameterSolution.compOrderParam(y));
                writer.flush();
                
                if(i < 1) {
                    for(int j = 0; j < y.length; ++j) {
                        y0[j] = y[j];
                    }
                }
            } catch(Exception ex) {
                writer.println(feff + ", 0");
                System.out.println(ex.getMessage());
            }
        }
        
        writer.close();

        long endTime = System.nanoTime();

        System.out.println("Run time: " + (endTime - startTime)/1.0e9 + " seconds");
    }
}
