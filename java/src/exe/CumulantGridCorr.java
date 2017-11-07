package exe;

import handlers.CumulantSteadyStateTerminator;
import handlers.WriteHandler;

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

public class CumulantGridCorr {

    /**
     * @param args
     * @throws FileNotFoundException 
     * @throws UnsupportedEncodingException 
     */
    public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException {
        
        String dir = "/Users/kristophertucker/output/discrete/wronggrid3/";
        int n = 70;
        double h = 0.0001;
        double gamma = 1.0;
        double tmax = 7.0;
        double tmin = 5.0;
        double f = 1.0;
        boolean correlate = true;
        
        // Get natural frequencies from Gaussian distribution
        double[] d = new double[n];
//        SynchUtils.detuneGauss(delta, d);
//        SynchUtils.detuneLor(delta, d);
        
        ArrayList<CumulantParams> params = new ArrayList<CumulantParams>();

        double dmin = 0.0;
        double dmax = 20.0;
        double dd = 1.0;
        
        double gmin = 0.0;
        double gmax = 10.0;
        double dg = 1.0;
        
        // Override default values with command line arguments
        if(args.length > 0) {
            if(args.length != 7) {
                System.out.println("Usage: CumulantGridCorr g0 dg gf d0 dd df outdir");
                return;
            }

            gmin = Double.parseDouble(args[0]);
            dg = Double.parseDouble(args[1]);
            gmax = Double.parseDouble(args[2]);
            
            dmin = Double.parseDouble(args[3]);
            dd = Double.parseDouble(args[4]);
            dmax = Double.parseDouble(args[5]);
            
            dir = args[6];
        }
        
//        Scanner wstream = new Scanner(new File("/Users/kristophertucker/Google Drive/Research/Synch/cumulant_all/wopt_twoatom_dlo.txt"));
//        wstream.useDelimiter("\n");
        for(double di = dmin; di <= dmax; di += dd) {
            SynchUtils.detuneDiscrete(di, d);
//            String[] line = wstream.next().split("\\s+");
//            int gidx = 1;
            for(double gi = gmin; gi <= gmax; gi += dg) {
                DynaComplex alpha = new DynaComplex(f, gi);
//                double w = Double.parseDouble(line[gidx]);
                double w = SynchUtils.getWOpt(n);
//                ++gidx;
//                System.out.print(w + " ");
                CumulantParams p = new CumulantParams(n, gamma, w, di, alpha, d);
                params.add(p);
            }
//            System.out.println("");
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

        File fdir = new File(dir);
        fdir.mkdirs();
        boolean success = true;
        for(int idx = 0; idx < params.size(); ++idx) {
            CumulantParams cparams = params.get(idx);
            CumulantAllToAllODEs codes = new CumulantAllToAllODEs(cparams);
            DynaComplexODEAdapter odes = new DynaComplexODEAdapter(codes);
            
            int[] out_col = {0, 1, 2, 3,
                    codes.getStartIdx(1), codes.getStartIdx(1) + 1, codes.getStartIdx(1) + 2, codes.getStartIdx(1) + 3, 
                    codes.getStartIdx(2), codes.getStartIdx(2) + 1, codes.getStartIdx(2) + 2, codes.getStartIdx(2) + 3,
                    codes.getStartIdx(3), codes.getStartIdx(3) + 1,
                    codes.getStartIdx(4), codes.getStartIdx(4) + 1,
                    codes.getStartIdx(5), codes.getStartIdx(5) + 1};
            
//            WriteHandlerCorr writeHandler = new WriteHandlerCorr(dir + "avg_" + cparams.getFilename(), n);
            WriteHandler writeHandler = new WriteHandler(dir + cparams.getFilename(), out_col, false);
            CumulantSteadyStateTerminator term = new CumulantSteadyStateTerminator(tmin, 0.015, 50, 1000000, 0.0025, cparams.getN());
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
                SynchUtils.compCorr(cparams, y, dir + "time_corr_" + cparams.getFilename(), SynchUtils.CorrelationType.THIRD);
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
