package exe;

import handlers.CumulantSteadyStateTerminator;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;

import ode.CumulantAllToAllODEs;
import ode.DynaComplexODEAdapter;

import org.apache.commons.math3.ode.nonstiff.AdamsMoultonIntegrator;

import utils.DynaComplex;
import utils.SynchUtils;

public class CumulantSearch {

    /**
     * @param args
     * @throws UnsupportedEncodingException 
     * @throws FileNotFoundException 
     */
    public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException {
        int n = 100;
        double h = 0.005;
        double w = 4.0;
        double gamma = 1.0;
        double tmax = 10.0;
        double delta = 7.5;
        
        double gmin = 0.0;
        double gmax = 10.0;
        double dg = 0.5;
        
        // Get natural frequencies from file
        double[] d = new double[n];
        SynchUtils.detuneFile("/Users/kristophertucker/Google Drive/Research/Synch/data/dels_all.txt", d);

        // Get natural frequencies from Gaussian distribution
//        SynchUtils.detuneGauss(delta, d);
        
        // Get natural frequencies from Cauchy distribution
//        TDistribution tdist = new TDistribution(new Well19937c(1), 1);
//        double[] d = tdist.sample(n);
//        for(int i = 0; i < n; ++i) {
//            d[i] = d[i]*delta;
//            //System.out.println(d[i]);
//        }
        
        // Set natural frequencies to zero
//        double[] d = new double[n];
//        for(int i = 0; i < n; ++i) {
//            d[i] = 1.0;
//        }
        
        int dim = SynchUtils.getDimension(n);
        
        DynaComplex[] z0 = new DynaComplex[dim];
//        for(int i = 0; i < dim; ++i) {
//            z0[i] = new DynaComplex(0, 0);
//        }
//        
//        for(int i = 0; i < n; ++i) {
//            z0[i] = new DynaComplex(0.2, 0.0);
//        }
        
        SynchUtils.initialize(z0, Math.PI/2.0, n);
        
        double[] y0 = new double[2*dim];
        DynaComplexODEAdapter.toReal(z0, y0);
        
//        WriteHandlerCorr writeHandler = new WriteHandlerCorr("/Users/kristophertucker/output/corr/corrN100_D7p5_g10p0_2.txt", n);
//        CumulantDataRecorder recorder = new CumulantDataRecorder(0.025, tmax, n);
//        AdamsMoultonIntegrator integrator = new AdamsMoultonIntegrator(2, h*1.0e-4, h, 1.0e-3, 1.0e-2);
//        integrator.addStepHandler(writeHandler);
//        integrator.addStepHandler(recorder);

        double[] y = new double[2*dim];

        long startTime = System.nanoTime();

        PrintWriter writer = new PrintWriter("/Users/kristophertucker/output/search_out.txt", "UTF-8");
        for(double g = gmin; g <= gmax; g += dg) {
            DynaComplex alpha = new DynaComplex(1, g);
            System.out.println("g: " + g);
            
            CumulantAllToAllODEs codes = new CumulantAllToAllODEs(n, gamma, w, alpha, d);
            DynaComplexODEAdapter odes = new DynaComplexODEAdapter(codes);
            
            CumulantSteadyStateTerminator term = new CumulantSteadyStateTerminator(1.0, 0.015, 50, 1000000, 0.001, n);
            AdamsMoultonIntegrator integrator = new AdamsMoultonIntegrator(2, h*1.0e-4, h, 1.0e-3, 1.0e-2);
//            integrator.addStepHandler(writeHandler);
            integrator.addStepHandler(term.getDetector());
            integrator.addEventHandler(term, Double.POSITIVE_INFINITY, 1.0e-12, 100);

            integrator.integrate(odes, 0, y0, tmax, y);
            
            if(!term.getSteadyStateReached()) {
                System.out.println("WARNING: Failed to reach steady state for g = " + g);
            }

            writer.print(g + ", " + SynchUtils.compCorr(y, n).getReal() + "\n");
        }
        writer.close();
        
        long endTime = System.nanoTime();

        System.out.println("Run time: " + (endTime - startTime)/1.0e9 + " seconds");
    }
}