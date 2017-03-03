package exe;

import handlers.CumulantSteadyStateTerminator;
import handlers.WriteHandlerCorr;

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
        int n = 150;
        double h = 0.005;
        double gamma = 1.0;
        double tmax = 15.0;
        double delta = 7.5;
        
        double gmin = 0.0;
        double gmax = 50.0;
        double dg = 1.0;
        
        double w = SynchUtils.getWOpt(n);

        System.out.println("w: " + w);
        
        // Get natural frequencies from file
        double[] d = new double[n];
//        SynchUtils.detuneFile("/Users/kristophertucker/Google Drive/Research/Synch/data/dels_all.txt", d);

        // Get natural frequencies from Gaussian distribution
        SynchUtils.detuneGauss(delta, d);
        
        int dim = SynchUtils.getDimension(n);
        
        DynaComplex[] z0 = new DynaComplex[dim];
        
        SynchUtils.initialize(z0, Math.PI/2.0, n);
        
        double[] y0 = new double[2*dim];
        DynaComplexODEAdapter.toReal(z0, y0);
        
        double[] y = new double[2*dim];

        long startTime = System.nanoTime();

        String dir = "/Users/kristophertucker/output/corr_opt/D7p5/run2/N" + n + "/"; 
        PrintWriter writer = new PrintWriter(dir + "search_out.txt", "UTF-8");
        for(double g = gmin; g <= gmax; g += dg) {
            DynaComplex alpha = new DynaComplex(1, g);
            System.out.println("g: " + g);
            String gStr = Double.toString(g).replace('.', 'p');
            
            CumulantAllToAllODEs codes = new CumulantAllToAllODEs(n, gamma, w, alpha, d);
            DynaComplexODEAdapter odes = new DynaComplexODEAdapter(codes);
            
            CumulantSteadyStateTerminator term = new CumulantSteadyStateTerminator(1.0, 0.015, 50, 1000000, 0.001, n);
            WriteHandlerCorr writeHandler = new WriteHandlerCorr(dir + "corrN" + n + "_g" + gStr + ".txt", n);
            AdamsMoultonIntegrator integrator = new AdamsMoultonIntegrator(2, h*1.0e-4, h, 1.0e-3, 1.0e-2);
            integrator.addStepHandler(writeHandler);
            integrator.addStepHandler(term.getDetector());
            integrator.addEventHandler(term, Double.POSITIVE_INFINITY, 1.0e-12, 100);

            integrator.integrate(odes, 0, y0, tmax, y);
            
            if(!term.getSteadyStateReached()) {
                System.out.println("WARNING: Failed to reach steady state for g = " + g);
                writer.print(g + ", " + -1.0 + ", " + term.getStopTime() + "\n");
            } else {
                writer.print(g + ", " + SynchUtils.compCorr(y, n).getReal() + ", " + term.getStopTime() + "\n");
                
            }

            writer.flush();
        }
        writer.close();
        
        long endTime = System.nanoTime();

        System.out.println("Run time: " + (endTime - startTime)/1.0e9 + " seconds");
    }
}