package exe;

import handlers.CumulantSteadyStateTerminator;
import handlers.WriteBlochVectors;
import handlers.WriteHandlerCorr;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.Scanner;

import ode.CumulantAllToAllODEs;
import ode.SystemParams;
import ode.DynaComplexODEAdapter;

import org.apache.commons.math3.ode.nonstiff.AdamsMoultonIntegrator;

import utils.DynaComplex;
import utils.SynchUtils;
import eval.CumulantEval;
import eval.SystemEval.InitAngleType;

public class CumulantSingle {

    /**
     * @param args
     * @throws UnsupportedEncodingException 
     * @throws FileNotFoundException 
     */
    public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException {
        int n = 2;
        double h = 0.001;
        double gamma = 1.0;
        double tmax = 0.25;
        double tmin = 0.25;
        double delta = 0.0;
        double f = 1.0;
        double g = 30.0;
        double gel = 0.0;
        boolean correlate = false;
        boolean outputBloch = false;
        boolean upper = false;

        double w = SynchUtils.getWOpt(n);
        w = 3.0;
        //w = 60.0;
        
        // Get natural frequencies from Gaussian distribution
        double[] d = new double[n];
        SynchUtils.detuneGauss(delta, d);
//        SynchUtils.detuneLor(delta, d);
//        SynchUtils.detuneDiscrete(delta, d);
        
        // Get natural frequencies from file
//        Scanner inputStream = new Scanner(new File("/Users/kristophertucker/Google Drive/Research/Synch/cumulant_all/even_delta_D100.txt"));
//        int idx = 0;
//        while(inputStream.hasNext()) {
//            d[idx] = Double.parseDouble(inputStream.next());
//            ++idx;
//        }
        
        DynaComplex alpha = new DynaComplex(f, g);
        
        double[] fvec = new double[n];
        double[] gvec = new double[n];
        for(int i = 0; i < n; ++i) {
            fvec[i] = f*Math.cos(7.0*Math.PI*i/6.0);
            gvec[i] = g*Math.cos(7.0*Math.PI*i/6.0);
        }
        
        SystemParams params = new SystemParams(n, gamma, w, delta, alpha, d);
        params.setGel(gel);
        
        CumulantAllToAllODEs codes = new CumulantAllToAllODEs(params);
//        CumulantDecoupledODEs codes = new CumulantDecoupledODEs(n, gamma, w, fvec, gvec, d);
        DynaComplexODEAdapter odes = new DynaComplexODEAdapter(codes);
        int dim = codes.getDimension();
        
        DynaComplex[] z0 = new DynaComplex[dim];
        SynchUtils.initialize(z0, Math.PI/2.0, n);
        
        // Get initial conditions from a file
        if(upper) {
            Scanner inputStream = new Scanner(new File("/Users/kristophertucker/output/grid/glow/upper/N30/D0p0/g10p0/final_w16p75.txt"));
            inputStream.useDelimiter("\n");
            int idx = 0;
            while(inputStream.hasNext()) {
                String[] line = inputStream.next().split(",");
                z0[idx] = new DynaComplex(Double.parseDouble(line[0]), Double.parseDouble(line[1]));
                if(z0[idx].mod() < 1.0e-10) {
                    z0[idx].set(0,0);
                }
                ++idx;
            }
        }
        
//        for(int i = 0; i < z0.length; ++i) {
//            System.out.println(z0[i]);
//        }
//        
//        if(n > 0) return;

        double[] y0 = new double[2*dim];
        DynaComplexODEAdapter.toReal(z0, y0);
        
        CumulantEval eval = new CumulantEval(n);
        eval.initialize(y0, 0.25*Math.PI, 0.25*Math.PI, InitAngleType.CONST, InitAngleType.CONST);
        
        int startIdx[] = new int[6];
        SynchUtils.getStartIdx(startIdx, n);
        double kickval = 1.0e-21;
        DynaComplex kick = new DynaComplex(0.182297,0.0334138);
        kick.multiply(kickval);
        for(int i = startIdx[0]; i < startIdx[1]; ++i) {
            //y0[2*i] += kickval;
        }
        
        for(int i = startIdx[1]; i < startIdx[2]; ++i) {
            //y0[2*i] += kick.getReal();
            //y0[2*i+1] += kick.getImaginary();
            //System.out.println(y0[2*i] + " " + y0[2*i+1]);
        }
        
        System.out.println("init corr: " + SynchUtils.compCorr(y0, n));

        double[] y = new double[2*dim];
        
        int[] out_col = {0, 1, 2, 3,
                codes.getStartIdx(1), codes.getStartIdx(1) + 1, codes.getStartIdx(1) + 2, codes.getStartIdx(1) + 3, 
                codes.getStartIdx(2), codes.getStartIdx(2) + 1, codes.getStartIdx(2) + 2, codes.getStartIdx(2) + 3,
                codes.getStartIdx(3), codes.getStartIdx(3) + 1,
                codes.getStartIdx(4), codes.getStartIdx(4) + 1,
                codes.getStartIdx(5), codes.getStartIdx(5) + 1};
        
        long startTime = System.nanoTime();

        String dir = "/Users/kristophertucker/output/temp/";
        if(upper) dir += "upper/";
        File fdir = new File(dir);
        fdir.mkdirs();
        WriteBlochVectors writeBloch = null;
        if(outputBloch) writeBloch = new WriteBlochVectors(dir + "bloch_" + params.getFilename(), n);
        WriteHandlerCorr writeHandler = new WriteHandlerCorr(dir + "avg_" + params.getFilename(), n);
//        WriteHandler writeHandler = new WriteHandler(dir + params.getFilename(), out_col, false);
        CumulantSteadyStateTerminator term = new CumulantSteadyStateTerminator(tmin, 0.015, 50, 1000000, 0.002, n);
        AdamsMoultonIntegrator integrator = new AdamsMoultonIntegrator(2, h*1.0e-4, h, 1.0e-3, 1.0e-2);
        //GraggBulirschStoerIntegrator integrator = new GraggBulirschStoerIntegrator(1.0e-18, h, 1.0e-3, 1.0e-2);
        //DormandPrince54Integrator integrator = new DormandPrince54Integrator(1.0e-18, h, 1.0e-3, 1.0e-2);
        //ClassicalRungeKuttaIntegrator integrator = new ClassicalRungeKuttaIntegrator(h);
        if(outputBloch) integrator.addStepHandler(writeBloch);
        integrator.addStepHandler(writeHandler);
        integrator.addStepHandler(term.getDetector());
        integrator.addEventHandler(term, Double.POSITIVE_INFINITY, 1.0e-12, 100);

        System.out.println("w: " + w);
        integrator.integrate(odes, 0, y0, tmax, y);
        
        System.out.println("final corr: " + SynchUtils.compCorr(y, n));
        
        if(correlate) {
            SynchUtils.compCorr(params, y, dir + "two_time_" + params.getFilename(), SynchUtils.CorrelationType.THIRD);
        }

        long endTime = System.nanoTime();
        
        DynaComplexODEAdapter.toComplex(y, z0);
        PrintWriter writer = new PrintWriter(dir + "final_answer_" + params.getFilename(), "UTF-8");
        for(int i = 0; i < dim; ++i) {
            writer.write(z0[i].getReal() + ", " + z0[i].getImaginary() + "\n");
        }
        writer.close();

        System.out.println("Run time: " + (endTime - startTime)/1.0e9 + " seconds");
    }
}
