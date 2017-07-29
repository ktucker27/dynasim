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

public class CumulantRun {

    /**
     * @param args
     * @throws FileNotFoundException 
     * @throws UnsupportedEncodingException 
     */
    public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException {
        
        int n = 120;
        double h = 0.0001;
        double gamma = 1.0;
        double tmax = 7.0;
        double tmin = 5.0;
        double delta = 0.0;
        double f = 1;
        double g = 3.0;
        double gel = 0.0;
        boolean correlate = false;
        
        int nmax = n;
//        int nmax = 30;
//        int nmin = 16;

        DynaComplex alpha = new DynaComplex(f, g);
        
        // Get natural frequencies from Gaussian distribution
        double[] d = new double[nmax];
//        SynchUtils.detuneGauss(delta, d);
//        SynchUtils.detuneLor(delta, d);
        SynchUtils.detuneDiscrete(delta, d);
        
        // Get natural frequencies from file
//        Scanner inputStream = new Scanner(new File("/Users/kristophertucker/Google Drive/Research/Synch/cumulant_all/even_delta_D100.txt"));
//        int didx = 0;
//        while(inputStream.hasNext()) {
//            d[didx] = Double.parseDouble(inputStream.next());
//            ++didx;
//            if(didx >= n) break;
//        }
//        
//        for(int i = 0; i < d.length; ++i) {
//            System.out.println(d[i]);
//        }
        
        // Get natural frequencies from Cauchy distribution
//        double delta = 2.0;
//        TDistribution tdist = new TDistribution(new Well19937c(1), 1);
//        double[] d = tdist.sample(n);
//        for(int i = 0; i < n; ++i) {
//            d[i] = d[i]*delta;
//            //System.out.println(d[i]);
//        }
        
        // Set natural frequencies to zero
//        double[] d = new double[n];
//        for(int i = 0; i < n; ++i) {
//            d[i] = 0.0;
//        }
        
        ArrayList<CumulantParams> params = new ArrayList<CumulantParams>();

        double wmin = 2.5;
        double wmax = 120.0;
        double dw = (wmax - wmin)/50;
        for(double w = wmin; w <= wmax; w += dw) {
            CumulantParams p = new CumulantParams(n, gamma, w, delta, alpha, d);
            p.setGel(gel);
            params.add(p);
        }
        
//        for(int nn = nmin; nn <= nmax; ++nn) {
//            CumulantParams p = new CumulantParams(nn, gamma, SynchUtils.getWOpt(nn), delta, alpha, d);
//            params.add(p);
//        }
        
        //DynaCumulantODEs codes = new DynaCumulantODEs(n, gamma, w, coupling, d);
//        CumulantAllToAllODEs codes = new CumulantAllToAllODEs(n, gamma, 0, alpha, d);
        int dim = SynchUtils.getDimension(nmax);
        
        DynaComplex[] z0 = new DynaComplex[dim];
        SynchUtils.initialize(z0, Math.PI/2.0, n);
        
        // Get initial conditions from a file
//        Scanner inputStream = new Scanner(new File("/Users/kristophertucker/output/twotime/long/backward/N30/D15p0/g20p0/final_w16p75.txt"));
//        inputStream.useDelimiter("\n");
//        int iidx = 0;
//        while(inputStream.hasNext()) {
//            String[] line = inputStream.next().split(",");
//            z0[iidx] = new DynaComplex(Double.parseDouble(line[0]), Double.parseDouble(line[1]));
//            if(z0[iidx].mod() < 1.0e-10) {
//                z0[iidx].set(0,0);
//            }
//            ++iidx;
//        }

        double[] y0 = new double[2*dim];
        DynaComplexODEAdapter.toReal(z0, y0);
        
//        double[] y00 = new double[2*dim];
//        DynaComplexODEAdapter.toReal(z0, y00);
//
//        double[] y0 = new double[2*SynchUtils.getDimension(params.get(0).getN())];
//        SynchUtils.changeDim(y00, y0, nmax, params.get(0).getN());

//        int n2 = n + 1;
//        double[] y1 = new double[2*SynchUtils.getDimension(n2)];
//        SynchUtils.changeDim(y0, y1, n, n2);
//        DynaComplex[] z1 = new DynaComplex[SynchUtils.getDimension(n2)];
//        for(int i = 0; i < z1.length; ++i) {
//            z1[i] = new DynaComplex();
//        }
//        DynaComplexODEAdapter.toComplex(y1, z1);
//        for(int i = 0; i < z1.length; ++i) {
//            System.out.println(z1[i]);
//        }
//        if(n > 0) return;

        System.out.println(params.get(0).toString());
        
//        int[] out_col = {0, 1, 2, 3,
//                        codes.getStartIdx(1), codes.getStartIdx(1) + 1, codes.getStartIdx(1) + 2, codes.getStartIdx(1) + 3, 
//                        codes.getStartIdx(2), codes.getStartIdx(2) + 1, codes.getStartIdx(2) + 2, codes.getStartIdx(2) + 3,
//                        codes.getStartIdx(3), codes.getStartIdx(3) + 1,
//                        codes.getStartIdx(4), codes.getStartIdx(4) + 1,
//                        codes.getStartIdx(5), codes.getStartIdx(5) + 1};
        
        long startTime = System.nanoTime();

        String dir = "/Users/kristophertucker/output/vw/" + params.get(0).getResultsDir().getAbsolutePath() + "/";
        File fdir = new File(dir);
        fdir.mkdirs();
        PrintWriter corrWriter = new PrintWriter(dir + "corr.txt", "UTF-8");
        boolean success = true;
        for(int idx = 0; idx < params.size(); ++idx) {
            CumulantParams cparams = params.get(idx);
            CumulantAllToAllODEs codes = new CumulantAllToAllODEs(cparams);
            DynaComplexODEAdapter odes = new DynaComplexODEAdapter(codes);

            String wStr = Double.toString(cparams.getW()).replace('.', 'p');
            
            WriteHandlerCorr writeHandler = new WriteHandlerCorr(dir + "avg_w" + wStr + ".txt", n);
            writeHandler.setMinTime(tmin - 5.0);
            CumulantSteadyStateTerminator term = new CumulantSteadyStateTerminator(tmin, 0.015, 50, 1000000, 0.0025, cparams.getN());
            AdamsMoultonIntegrator integrator = new AdamsMoultonIntegrator(2, h*1.0e-4, h, 1.0e-3, 1.0e-2);
            //GraggBulirschStoerIntegrator integrator = new GraggBulirschStoerIntegrator(1.0e-18, h, 1.0e-3, 1.0e-2);
            //DormandPrince54Integrator integrator = new DormandPrince54Integrator(1.0e-18, h, 1.0e-3, 1.0e-2);
            //ClassicalRungeKuttaIntegrator integrator = new ClassicalRungeKuttaIntegrator(h);
            integrator.addStepHandler(writeHandler);
            integrator.addStepHandler(term.getDetector());
            integrator.addEventHandler(term, Double.POSITIVE_INFINITY, 1.0e-12, 100);

            double[] y = new double[2*SynchUtils.getDimension(cparams.getN())];
            System.out.println("w: " + cparams.getW());
            System.out.println("init corr: " + SynchUtils.compCorr(y0, cparams.getN()));
            integrator.integrate(odes, 0, y0, tmax, y);
            System.out.println("final corr: " + SynchUtils.compCorr(y, cparams.getN()));
            
            // Copy solution to initial conditions
//            for(int i = 0; i < y.length; ++i) {
//                if(Math.abs(y[i]) < 1.0e-10) {
//                    y0[i] = 0.0;
//                } else {
//                    y0[i] = y[i];
//                }
//            }
            
//            if(idx+1 < params.size()) {
//                y0 = new double[2*SynchUtils.getDimension(params.get(idx+1).getN())];
//                SynchUtils.changeDim(y00, y0, nmax, params.get(idx+1).getN());
//            }

            // Compute the correlation function if requested
            if(correlate) {
                SynchUtils.compCorr(cparams, y, dir + "time_corr_" + wStr + ".txt");
            }

            DynaComplexODEAdapter.toComplex(y, z0);
            PrintWriter writer = new PrintWriter(dir + "final_w" + wStr + ".txt", "UTF-8");
            for(int i = 0; i < SynchUtils.getDimension(cparams.getN()); ++i) {
                writer.write(z0[i].getReal() + ", " + z0[i].getImaginary() + "\n");
            }
            writer.close();
            
            if(!term.getSteadyStateReached()) {
                System.out.println("WARNING: Failed to reach steady state for w = " + cparams.getW());
                corrWriter.print(cparams.getW() + ", " + -1.0 + ", " + term.getStopTime() + "\n");
                success = false;
            } else {
                corrWriter.print(cparams.getW() + ", " + SynchUtils.compCorr(y, cparams.getN()).getReal() + ", " + term.getStopTime() + "\n");
            }
            corrWriter.flush();
        }
        corrWriter.close();
        
        long endTime = System.nanoTime();

        System.out.println(success);
        System.out.println("Run time: " + (endTime - startTime)/1.0e9 + " seconds");        
    }
}
