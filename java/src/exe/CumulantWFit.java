package exe;

import java.io.FileNotFoundException;
import java.io.UnsupportedEncodingException;

import ode.CumulantAllToAllODEs;
import ode.DynaComplexODEAdapter;

import org.apache.commons.math3.ode.nonstiff.AdamsMoultonIntegrator;

import utils.DynaComplex;
import utils.SynchUtils;
import utils.WriteHandlerCorr;

public class CumulantWFit {

    /**
     * @param args
     * @throws UnsupportedEncodingException 
     * @throws FileNotFoundException 
     */
    public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException {
        int n = 10;
        double h = 0.01;
        double wmin = 1.0;
        double wmax = 10.0;
        double gamma = 1.0;
        double tmax = 3.0;
        double delta = 7.5;
        
        double target = 0.01224;
        double tol = 1.0e-5;
        
        // Get natural frequencies from Gaussian distribution
        double[] d = new double[n];
        SynchUtils.detuneGauss(delta, d);
        
        DynaComplex alpha = new DynaComplex(1, 0.0);
        
        CumulantAllToAllODEs codes = new CumulantAllToAllODEs(n, gamma, 0, alpha, d);
        DynaComplexODEAdapter odes = new DynaComplexODEAdapter(codes);
        int dim = codes.getDimension();

        DynaComplex[] z0 = new DynaComplex[dim];
        
        SynchUtils.initialize(z0, Math.PI/2.0, n);
        
        double[] y0 = new double[2*dim];
        DynaComplexODEAdapter.toReal(z0, y0);
                
        double[] y = new double[2*dim];

        long startTime = System.nanoTime();

        double corr = 0.0;
        double minCorr = 0.0;
        double maxCorr = 0.0;
        double w = 0.0;
        boolean initMin = false;
        boolean initMax = false;
        while(!initMin || !initMax || Math.abs(corr - target) > tol) {
            System.out.println("w range: " + wmin + ", " + wmax);

            WriteHandlerCorr writeHandler = new WriteHandlerCorr("/Users/kristophertucker/output/corr/fit.txt", n);
            AdamsMoultonIntegrator integrator = new AdamsMoultonIntegrator(2, h*1.0e-4, h, 1.0e-3, 1.0e-2);
            integrator.addStepHandler(writeHandler);
          
//            SpinCorrSteadyStateTest test = new SpinCorrSteadyStateTest(1.0e-5, n);
//            SynchSteadyStateTerminator term = new SynchSteadyStateTerminator(test, 3.0, 200, 5000000);
//            AdamsMoultonIntegrator integrator = new AdamsMoultonIntegrator(2, h*1.0e-4, h, 1.0e-3, 1.0e-2);
//            integrator.addStepHandler(term.getDetector());
//            integrator.addEventHandler(term, Double.POSITIVE_INFINITY, 1.0e-12, 100);

            if(!initMin) {
                w = wmin;
            } else if(!initMax) {
                w = wmax;
            } else {
                w = 0.5*(wmin + wmax);
            }
            
            codes.setW(w);
            integrator.integrate(odes, 0, y0, tmax, y);
            
//            if(!term.getReturnOk()) {
//                System.out.println("Failed to find steady state for w = " + w);
//                break;
//            }
            
            corr = SynchUtils.compCorr(y, n).getReal();
            System.out.println("w = " + w + ": " + corr);

            if(!initMin) {
                initMin = true;
                minCorr = corr;
            } else if(!initMax) {
                initMax = true;
                maxCorr = corr;
                if(minCorr > target || maxCorr < target) {
                    System.out.println("Failed to bracket target with initial values");
                    break;
                }
            } else {
                if(corr < target) {
                    wmin = w;
                    minCorr = corr;
                } else {
                    wmax = w;
                    maxCorr = corr;
                }
            }
        }
        
        long endTime = System.nanoTime();

        System.out.println("Final w: " + w);
        System.out.println("Run time: " + (endTime - startTime)/1.0e9 + " seconds");
    }

}
