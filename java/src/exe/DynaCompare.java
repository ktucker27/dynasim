package exe;

import handlers.CumulantSteadyStateTerminator;

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

public class DynaCompare {

    /**
     * @param args
     * @throws FileNotFoundException 
     * @throws UnsupportedEncodingException 
     */
    public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException {
        
        int n = 16;
        double h = 0.005;
        double gamma = 1.0;
        double tmax = 20.0;
        double delta = 0.0;
        double f = 1.0;
        double g = 0.0;
        //boolean correlate = false;
        
        DynaComplex alpha = new DynaComplex(f, g);
        
        double wmin = 2.5;
        double wmax = 40.0;
        double dw = (wmax - wmin)/100;

        // Get natural frequencies from Gaussian distribution
        double[] d = new double[n];
        SynchUtils.detuneGauss(delta, d);
        
        // Get natural frequencies from file
//        double[] d = new double[n];
//        Scanner inputStream = new Scanner(new File("/Users/kristophertucker/Google Drive/Research/Synch/data/dels.txt"));
//        int idx = 0;
//        while(inputStream.hasNext()) {
//            d[idx] = Double.parseDouble(inputStream.next());
//            ++idx;
//            if(idx >= n) break;
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
        for(double w = wmin; w <= wmax; w += dw) {
            CumulantParams p = new CumulantParams(n, gamma, w, delta, alpha, d);
            params.add(p);
        }
        
        //DynaCumulantODEs codes = new DynaCumulantODEs(n, gamma, w, coupling, d);
        CumulantAllToAllODEs codes = new CumulantAllToAllODEs(n, gamma, 0, alpha, d);
        int dim = codes.getDimension();
        
        DynaComplex[] z0 = new DynaComplex[dim];
        SynchUtils.initialize(z0, Math.PI/2.0, n);
        
        // Get initial conditions from a file
//        Scanner inputStream = new Scanner(new File("/Users/kristophertucker/output/info/wruns_chao/g20/info_out_w120p0.txt"));
//        inputStream.useDelimiter("\n");
//        int idx = 0;
//        while(inputStream.hasNext()) {
//            String[] line = inputStream.next().split(",");
//            z0[idx] = new DynaComplex(Double.parseDouble(line[0]), Double.parseDouble(line[1]));
//            ++idx;
//        }

//        for(int i = 0; i < dim; ++i) {
//            z0[i] = new DynaComplex(0, 0);
//        }
        
//        for(int i = 0; i < n; ++i) {
//            z0[i] = new DynaComplex(0.2, 0.0);
//        }
        
        double[] y0 = new double[2*dim];
        DynaComplexODEAdapter.toReal(z0, y0);
        
        System.out.println("init corr: " + SynchUtils.compCorr(y0, n));
        
//        int[] out_col = {0, 1, 2, 3,
//                        codes.getStartIdx(1), codes.getStartIdx(1) + 1, codes.getStartIdx(1) + 2, codes.getStartIdx(1) + 3, 
//                        codes.getStartIdx(2), codes.getStartIdx(2) + 1, codes.getStartIdx(2) + 2, codes.getStartIdx(2) + 3,
//                        codes.getStartIdx(3), codes.getStartIdx(3) + 1,
//                        codes.getStartIdx(4), codes.getStartIdx(4) + 1,
//                        codes.getStartIdx(5), codes.getStartIdx(5) + 1};

        double[] y = new double[2*dim];
        
        long startTime = System.nanoTime();

        String dir = "/Users/kristophertucker/output/vw/" + params.get(0).getResultsDir().getAbsolutePath() + "/";
        File fdir = new File(dir);
        fdir.mkdirs();
        PrintWriter corrWriter = new PrintWriter(dir + "corr.txt", "UTF-8");
        boolean success = true;
        for(int idx = 0; idx < params.size(); ++idx) {
            CumulantParams cparams = params.get(idx);
            codes = new CumulantAllToAllODEs(cparams);
            DynaComplexODEAdapter odes = new DynaComplexODEAdapter(codes);
            
            //WriteHandlerCorr writeHandler = new WriteHandlerCorr(dir + "full.txt", n);
            CumulantSteadyStateTerminator term = new CumulantSteadyStateTerminator(1.0, 0.015, 50, 1000000, 0.0025, n);
            AdamsMoultonIntegrator integrator = new AdamsMoultonIntegrator(2, h*1.0e-4, h, 1.0e-3, 1.0e-2);
            //GraggBulirschStoerIntegrator integrator = new GraggBulirschStoerIntegrator(1.0e-18, h, 1.0e-3, 1.0e-2);
            //DormandPrince54Integrator integrator = new DormandPrince54Integrator(1.0e-18, h, 1.0e-3, 1.0e-2);
            //ClassicalRungeKuttaIntegrator integrator = new ClassicalRungeKuttaIntegrator(h);
            //integrator.addStepHandler(writeHandler);
            integrator.addStepHandler(term.getDetector());
            integrator.addEventHandler(term, Double.POSITIVE_INFINITY, 1.0e-12, 100);

            System.out.println("w: " + cparams.getW());
            integrator.integrate(odes, 0, y0, tmax, y);
            
            // Copy solution to initial conditions
//            for(int i = 0; i < y.length; ++i) {
//                y0[i] = y[i];
//            }

            DynaComplexODEAdapter.toComplex(y, z0);
            String wStr = Double.toString(cparams.getW()).replace('.', 'p');
            PrintWriter writer = new PrintWriter(dir + "final_w" + wStr + ".txt", "UTF-8");
            for(int i = 0; i < dim; ++i) {
                writer.write(z0[i].getReal() + ", " + z0[i].getImaginary() + "\n");
            }
            writer.close();
            
            if(!term.getSteadyStateReached()) {
                System.out.println("WARNING: Failed to reach steady state for w = " + cparams.getW());
                corrWriter.print(cparams.getW() + ", " + -1.0 + ", " + term.getStopTime() + "\n");
                success = false;
            } else {
                corrWriter.print(cparams.getW() + ", " + SynchUtils.compCorr(y, n).getReal() + ", " + term.getStopTime() + "\n");
            }
            corrWriter.flush();
        }
        corrWriter.close();
        
        long endTime = System.nanoTime();

        System.out.println(success);
        System.out.println("Run time: " + (endTime - startTime)/1.0e9 + " seconds");
        
        // Compute the correlation function if requested
        /*
        if(correlate) {
            DynaComplex[] z = new DynaComplex[dim];
            for(int i = 0; i < dim; ++i) {
                z[i] = new DynaComplex(0, 0);
            }
            DynaComplexODEAdapter.toComplex(y, z);
            
            DynaComplex[] szs = new DynaComplex[n];
            int idx2 = 0;
            for(int i = codes.getStartIdx(2)/2; i < codes.getStartIdx(2)/2 + n; ++i) {
                szs[idx2] = new DynaComplex(z[i].getReal(), z[i].getImaginary());
                ++idx2;
            }
            
            double mod_sum = 0.0;
            DynaComplex[] z02 = new DynaComplex[n*n];
            idx2 = codes.getStartIdx(3)/2;
            for(int i = 0; i < n; ++i) {
                z02[i*n + i] = new DynaComplex(1.0, 0);
                for(int j = i+1; j < n; ++j) {
                    z02[i*n + j] = new DynaComplex(z[idx2].getReal(), z[idx2].getImaginary());
                    z02[j*n + i] = new DynaComplex(z[idx2].getReal(), z[idx2].getImaginary());
                    z02[j*n + i].conjugate();
                    
                    mod_sum += 2.0*z02[j*n + i].mod();
                    
                    //System.out.println(z02[i*n + j].getReal() + "+j*" + z02[i*n + j].getImaginary() + " " + idx2 + " " + codes.getStartIdx(3)/2 + " " + codes.getStartIdx(4)/2);
                    ++idx2;
                }
            }
            System.out.println("Avg corr factor: " + mod_sum/(n*n));
            
            CorrelationODEs c_corr_odes = new CorrelationODEs(n, gamma, w, coupling, d, szs);
            odes = new DynaComplexODEAdapter(c_corr_odes);
            
            int[] out_col2 = new int[2*n*n];
            for(int i = 0; i < 2*n*n; ++i) {
                out_col2[i] = i;
            }
            
            //writeHandler = new WriteHandler("/Users/kristophertucker/output/correlation/corr_N32_al2p0_D7p5.txt", out_col2);
            AdamsMoultonIntegrator integrator = new AdamsMoultonIntegrator(2, 1.0e-18, .001, 1.0e-3, 1.0e-2);
            integrator.clearStepHandlers();
            //integrator.addStepHandler(writeHandler);
            
            double[] y02 = new double[2*n*n];
            DynaComplexODEAdapter.toReal(z02, y02);
            
            double[] y2 = new double[2*n*n];
            
            startTime = System.nanoTime();
            
            integrator.integrate(odes, 0, y02, 5, y2);
            
            endTime = System.nanoTime();
            
            System.out.println("Correlation time: " + (endTime - startTime)/1.0e9 + " seconds");
        }
        */
    }

}
