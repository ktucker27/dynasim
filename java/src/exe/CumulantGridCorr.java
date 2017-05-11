package exe;

import handlers.CumulantSteadyStateTerminator;
import handlers.WriteHandlerCorr;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Scanner;

import ode.CumulantAllToAllODEs;
import ode.CumulantParams;
import ode.DynaComplexODEAdapter;

import org.apache.commons.math3.ode.nonstiff.AdamsMoultonIntegrator;

import utils.DynaComplex;
import utils.SynchUtils;

public class CumulantGridCorr {

    /**
     * @param args
     * @throws FileNotFoundException 
     * @throws UnsupportedEncodingException 
     */
    public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException {
        
        int n = 30;
        double h = 0.001;
        double gamma = 1.0;
        double tmax = 20.0;
        double f = 1.0;
        boolean correlate = true;
        
        // Get natural frequencies from Gaussian distribution
        double[] d = new double[n];
//        SynchUtils.detuneGauss(delta, d);
//        SynchUtils.detuneLor(delta, d);
        
        ArrayList<CumulantParams> params = new ArrayList<CumulantParams>();

        double dmin = 2.0;
        double dmax = 20.0;
        double dd = 2.0;
        
        double gmin = 5.0;
        double gmax = 40.0;
        double dg = 5.0;
        
        Scanner wstream = new Scanner(new File("/Users/kristophertucker/Google Drive/Research/Synch/cumulant_all/wopt_grid.txt"));
        wstream.useDelimiter("\n");
        for(double di = dmin; di <= dmax; di += dd) {
            SynchUtils.detuneGauss(di, d);
            //String[] line = wstream.next().split("\\s+");
            int gidx = 1;
            for(double gi = gmin; gi <= gmax; gi += dg) {
                DynaComplex alpha = new DynaComplex(f, gi);
                //double w = Double.parseDouble(line[gidx]);
                double w = 16.75;
                ++gidx;
                System.out.print(w + " ");
                CumulantParams p = new CumulantParams(n, gamma, w, di, alpha, d);
                params.add(p);
            }
            System.out.println("");
        }
        
        int dim = SynchUtils.getDimension(n);
        
        DynaComplex[] z0 = new DynaComplex[dim];
        
        SynchUtils.initialize(z0, Math.PI/2.0, n);
        
        double[] y0 = new double[2*dim];
        DynaComplexODEAdapter.toReal(z0, y0);
        
        // Remove small y0 values to avoid underflow issues
        for(int i = 0; i < y0.length; ++i) {
            if(Math.abs(y0[i]) < 1.0e-10) {
                y0[i] = 0.0;
            }
        }
        
        long startTime = System.nanoTime();

        String dir = "/Users/kristophertucker/output/grid/twotime/";
        File fdir = new File(dir);
        fdir.mkdirs();
        boolean success = true;
        for(int idx = 0; idx < params.size(); ++idx) {
            CumulantParams cparams = params.get(idx);
            CumulantAllToAllODEs codes = new CumulantAllToAllODEs(cparams);
            DynaComplexODEAdapter odes = new DynaComplexODEAdapter(codes);
            
            WriteHandlerCorr writeHandler = new WriteHandlerCorr(dir + "corr_" + cparams.getFilename(), n);
            CumulantSteadyStateTerminator term = new CumulantSteadyStateTerminator(5.0, 0.015, 50, 1000000, 0.0025, cparams.getN());
            AdamsMoultonIntegrator integrator = new AdamsMoultonIntegrator(2, h*1.0e-4, h, 1.0e-3, 1.0e-2);
            integrator.addStepHandler(writeHandler);
            integrator.addStepHandler(term.getDetector());
            integrator.addEventHandler(term, Double.POSITIVE_INFINITY, 1.0e-12, 100);

            double[] y = new double[2*SynchUtils.getDimension(cparams.getN())];
            System.out.print(cparams.toString());
            System.out.println("init corr: " + SynchUtils.compCorr(y0, cparams.getN()));
            integrator.integrate(odes, 0, y0, tmax, y);
            System.out.println("final corr: " + SynchUtils.compCorr(y, cparams.getN()) + "\n");
            
            // Compute the correlation function if requested
            if(correlate) {
                SynchUtils.compCorr(cparams, y, dir + "time_corr_" + cparams.getFilename());
            }

            DynaComplexODEAdapter.toComplex(y, z0);
            PrintWriter writer = new PrintWriter(dir + "final_w_" + cparams.getFilename(), "UTF-8");
            for(int i = 0; i < SynchUtils.getDimension(cparams.getN()); ++i) {
                writer.write(z0[i].getReal() + ", " + z0[i].getImaginary() + "\n");
            }
            writer.close();
            
            if(!term.getSteadyStateReached()) {
                System.out.println("WARNING: Failed to reach steady state for w = " + cparams.getW());
                success = false;
            }
        }
        
        long endTime = System.nanoTime();

        System.out.println(success);
        System.out.println("Run time: " + (endTime - startTime)/1.0e9 + " seconds");        
    }
}
