package exe;

import handlers.CumulantSteadyStateTerminator;
import handlers.WriteHandlerCorr;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;

import ode.CumulantAllToAllODEs;
import ode.CumulantParams;
import ode.DynaComplexODEAdapter;

import org.apache.commons.math3.ode.nonstiff.AdamsMoultonIntegrator;

import utils.DynaComplex;
import utils.SynchUtils;

public class CumulantGrid {

    /**
     * @param args
     * @throws FileNotFoundException 
     * @throws UnsupportedEncodingException 
     */
    public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException {
        
        int n = 70;
        double h = 0.0001;
        double gamma = 1.0;
        double tmax = 7.0;
        double tmin = 5.0;
        double f = 1.0;
        boolean correlate = false;
        boolean upper = false;
        
        // Get natural frequencies from Gaussian distribution
        double[] d = new double[n];
//        SynchUtils.detuneGauss(delta, d);
//        SynchUtils.detuneLor(delta, d);
        
        ArrayList<CumulantParams> params = new ArrayList<CumulantParams>();

        double wmin = 2.5;
        double wmax = 120.0;
        double dw = (wmax - wmin)/50;
        
        double dmin = 0.0;
        double dmax = 0.0;
        double dd = 2.0;
        
        double gmin = 20.0;
        double gmax = 20.0;
        double dg = 2.0;
                
        for(double di = dmin; di <= dmax; di += dd) {
//            SynchUtils.detuneGauss(di, d);
            SynchUtils.detuneDiscrete(di, d);
            for(double gi = gmin; gi <= gmax; gi += dg) {
                DynaComplex alpha = new DynaComplex(f, gi);
                for(double w = wmax; w >= wmin; w -= dw) {
                    CumulantParams p = new CumulantParams(n, gamma, w, di, alpha, d);
                    params.add(p);
                }
            }
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

        String dir;
        if(upper) {
            dir = "/Users/kristophertucker/output/discrete/grid2/upper/";
        } else {
            dir = "/Users/kristophertucker/output/discrete/ctest";
        }

        String prevdir = "";
        PrintWriter corrWriter = null;
        boolean success = true;
        for(int idx = 0; idx < params.size(); ++idx) {
            CumulantParams cparams = params.get(idx);
            CumulantAllToAllODEs codes = new CumulantAllToAllODEs(cparams);
            DynaComplexODEAdapter odes = new DynaComplexODEAdapter(codes);
            
            String resdir = dir + cparams.getResultsDir().getAbsolutePath() + "/";
            if(!resdir.equals(prevdir)) {
                System.out.println("Creating directory " + resdir);
                File fdir = new File(resdir);
                fdir.mkdirs();
                corrWriter = new PrintWriter(resdir + "corr.txt", "UTF-8");
                prevdir = resdir;
            }
            
            String wStr = Double.toString(cparams.getW()).replace('.', 'p');
            WriteHandlerCorr writeHandler = new WriteHandlerCorr(resdir + "avg_w" + wStr + ".txt", n);
            writeHandler.setMinTime(tmin - 2.0);
            //WriteHandlerCorr writeHandler = new WriteHandlerCorr(dir + "full.txt", n);
            CumulantSteadyStateTerminator term = new CumulantSteadyStateTerminator(tmin, 0.015, 50, 1000000, 0.0025, cparams.getN());
            AdamsMoultonIntegrator integrator = new AdamsMoultonIntegrator(2, h*1.0e-4, h, 1.0e-3, 1.0e-2);
            integrator.addStepHandler(writeHandler);
            integrator.addStepHandler(term.getDetector());
            integrator.addEventHandler(term, Double.POSITIVE_INFINITY, 1.0e-12, 100);

            double[] y = new double[2*SynchUtils.getDimension(cparams.getN())];
            System.out.println(cparams.toString());
            System.out.println("init corr: " + SynchUtils.compCorr(y0, cparams.getN()));
            integrator.integrate(odes, 0, y0, tmax, y);
            System.out.println("final corr: " + SynchUtils.compCorr(y, cparams.getN()));
            
            // Copy solution to initial conditions
            if(upper) {
                for(int i = 0; i < y.length; ++i) {
                    if(Math.abs(y[i]) < 1.0e-10) {
                        y0[i] = 0.0;
                    } else {
                        y0[i] = y[i];
                    }
                }
            }
            
            // Compute the correlation function if requested
            if(correlate) {
                SynchUtils.compCorr(cparams, y, resdir + "time_corr_" + wStr + ".txt");
            }

            DynaComplexODEAdapter.toComplex(y, z0);
            PrintWriter writer = new PrintWriter(resdir + "final_w" + wStr + ".txt", "UTF-8");
            for(int i = 0; i < SynchUtils.getDimension(cparams.getN()); ++i) {
                writer.write(z0[i].getReal() + ", " + z0[i].getImaginary() + "\n");
            }
            writer.close();
            
            if(!term.getSteadyStateReached()) {
                System.out.println("WARNING: Failed to reach steady state for w = " + cparams.getW());
                corrWriter.print(cparams.getW() + ", " + -1.0 + ", " + term.getStopTime() + "\n");
                success = false;
            } else {
                DynaComplex sigmap = SynchUtils.compSigmapAvg(y, cparams.getN());
                corrWriter.print(cparams.getW() + ", " + SynchUtils.compCorr(y, cparams.getN()).getReal());
                corrWriter.print(", " + sigmap.getReal());
                corrWriter.print(", " + sigmap.getImaginary());
                corrWriter.print(", " + term.getStopTime() + "\n");
            }
            corrWriter.flush();
        }
        corrWriter.close();
        
        long endTime = System.nanoTime();

        System.out.println(success);
        System.out.println("Run time: " + (endTime - startTime)/1.0e9 + " seconds");        
    }
}
