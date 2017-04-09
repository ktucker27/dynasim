package exe;

import handlers.CumulantSteadyStateTerminator;
import handlers.WriteHandlerCorr;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.Scanner;

import ode.CumulantAllToAllODEs;
import ode.DynaComplexODEAdapter;

import org.apache.commons.math3.ode.nonstiff.AdamsMoultonIntegrator;

import utils.DynaComplex;
import utils.SynchUtils;

public class CumulantSingle {

    /**
     * @param args
     * @throws UnsupportedEncodingException 
     * @throws FileNotFoundException 
     */
    public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException {
        int n = 16;
        double h = 0.001;
        double gamma = 1.0;
        double tmax = 10.0;
        double delta = 15.0;
        double f = 1.0;
        double g = 20.0;

        double w = SynchUtils.getW_D0(n);
        w = 10.0;
        
        // Get natural frequencies from Gaussian distribution
        double[] d = new double[n];
        SynchUtils.detuneGauss(delta, d);
        
        // Get natural frequencies from file
//        Scanner inputStream = new Scanner(new File("/Users/kristophertucker/Google Drive/Research/Synch/cumulant_all/even_delta_D100.txt"));
//        int idx = 0;
//        while(inputStream.hasNext()) {
//            d[idx] = Double.parseDouble(inputStream.next());
//            ++idx;
//        }
        
        for(int i = 0; i < d.length; ++i) {
            System.out.println(d[i]);
        }
        
        if(n > 0) {
            return;
        }
        
        DynaComplex alpha = new DynaComplex(f, g);
        
        double[] fvec = new double[n];
        double[] gvec = new double[n];
        for(int i = 0; i < n; ++i) {
            fvec[i] = f*Math.cos(7.0*Math.PI*i/6.0);
            gvec[i] = g*Math.cos(7.0*Math.PI*i/6.0);
        }
        
        CumulantAllToAllODEs codes = new CumulantAllToAllODEs(n, gamma, w, alpha, d);
//        CumulantDecoupledODEs codes = new CumulantDecoupledODEs(n, gamma, w, fvec, gvec, d);
        DynaComplexODEAdapter odes = new DynaComplexODEAdapter(codes);
        int dim = codes.getDimension();
        
        DynaComplex[] z0 = new DynaComplex[dim];
        SynchUtils.initialize(z0, Math.PI/2.0, n);

        // Get initial conditions from a file
//        Scanner inputStream = new Scanner(new File("/Users/kristophertucker/output/vw/long/forward/N30/D15p0/g20p0/final_w13p0.txt"));
//        inputStream.useDelimiter("\n");
//        int idx = 0;
//        while(inputStream.hasNext()) {
//            String[] line = inputStream.next().split(",");
//            z0[idx] = new DynaComplex(Double.parseDouble(line[0]), Double.parseDouble(line[1]));
//            ++idx;
//        }
        
        double[] y0 = new double[2*dim];
        DynaComplexODEAdapter.toReal(z0, y0);
        
        System.out.println("init corr: " + SynchUtils.compCorr(y0, n));

        double[] y = new double[2*dim];
        
        long startTime = System.nanoTime();

        String dir = "/Users/kristophertucker/output/temp/";
        WriteHandlerCorr writeHandler = new WriteHandlerCorr(dir + "corr_v_time.txt", n);
        CumulantSteadyStateTerminator term = new CumulantSteadyStateTerminator(1.0, 0.015, 50, 1000000, 0.002, n);
        AdamsMoultonIntegrator integrator = new AdamsMoultonIntegrator(2, h*1.0e-4, h, 1.0e-3, 1.0e-2);
        //GraggBulirschStoerIntegrator integrator = new GraggBulirschStoerIntegrator(1.0e-18, h, 1.0e-3, 1.0e-2);
        //DormandPrince54Integrator integrator = new DormandPrince54Integrator(1.0e-18, h, 1.0e-3, 1.0e-2);
        //ClassicalRungeKuttaIntegrator integrator = new ClassicalRungeKuttaIntegrator(h);
        integrator.addStepHandler(writeHandler);
        integrator.addStepHandler(term.getDetector());
        integrator.addEventHandler(term, Double.POSITIVE_INFINITY, 1.0e-12, 100);

        System.out.println("w: " + w);
        integrator.integrate(odes, 0, y0, tmax, y);

        long endTime = System.nanoTime();
        
        DynaComplexODEAdapter.toComplex(y, z0);
        PrintWriter writer = new PrintWriter(dir + "final_answer.txt", "UTF-8");
        for(int i = 0; i < dim; ++i) {
            writer.write(z0[i].getReal() + ", " + z0[i].getImaginary() + "\n");
        }
        writer.close();

        System.out.println("Run time: " + (endTime - startTime)/1.0e9 + " seconds");
    }

}
