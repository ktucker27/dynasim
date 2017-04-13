package exe;

import handlers.CumulantSteadyStateTerminator;
import handlers.TwoTimeHandler;
import handlers.WriteHandlerCorr;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;

import ode.CorrelationODEs;
import ode.CumulantAllToAllODEs;
import ode.CumulantParams;
import ode.DynaComplexODEAdapter;

import org.apache.commons.math3.ode.nonstiff.AdamsMoultonIntegrator;

import utils.DynaComplex;
import utils.SynchUtils;
import coupling.DynaConstCoupling;

public class CumulantSingle {

    /**
     * @param args
     * @throws UnsupportedEncodingException 
     * @throws FileNotFoundException 
     */
    public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException {
        int n = 30;
        double h = 0.001;
        double gamma = 1.0;
        double tmax = 10.0;
        double delta = 16.0;
        double f = 1.0;
        double g = 100.0;
        boolean correlate = true;

        double w = SynchUtils.getWOpt(n);
        
        // Get natural frequencies from Gaussian distribution
        double[] d = new double[n];
//        SynchUtils.detuneGauss(delta, d);
        SynchUtils.detuneLor(delta, d);
        
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
        
        CumulantParams params = new CumulantParams(n, gamma, w, delta, alpha, d);
        
        CumulantAllToAllODEs codes = new CumulantAllToAllODEs(params);
//        CumulantDecoupledODEs codes = new CumulantDecoupledODEs(n, gamma, w, fvec, gvec, d);
        DynaComplexODEAdapter odes = new DynaComplexODEAdapter(codes);
        int dim = codes.getDimension();
        
        DynaComplex[] z0 = new DynaComplex[dim];
        SynchUtils.initialize(z0, Math.PI/2.0, n);
        
        // Get initial conditions from a file
//        Scanner inputStream = new Scanner(new File("/Users/kristophertucker/output/twotime/long/backward/N30/D15p0/g20p0/final_w16p75.txt"));
//        inputStream.useDelimiter("\n");
//        int idx = 0;
//        while(inputStream.hasNext()) {
//            String[] line = inputStream.next().split(",");
//            z0[idx] = new DynaComplex(Double.parseDouble(line[0]), Double.parseDouble(line[1]));
//            if(z0[idx].mod() < 1.0e-10) {
//                z0[idx].set(0,0);
//            }
//            ++idx;
//        }
        
//        for(int i = 0; i < z0.length; ++i) {
//            System.out.println(z0[i]);
//        }
//        
//        if(n > 0) return;

        double[] y0 = new double[2*dim];
        DynaComplexODEAdapter.toReal(z0, y0);
        
        System.out.println("init corr: " + SynchUtils.compCorr(y0, n));

        double[] y = new double[2*dim];
        
        long startTime = System.nanoTime();

        String dir = "/Users/kristophertucker/output/temp/";
        WriteHandlerCorr writeHandler = new WriteHandlerCorr(dir + "corr_v_time.txt", n);
        CumulantSteadyStateTerminator term = new CumulantSteadyStateTerminator(2.0, 0.015, 50, 1000000, 0.002, n);
        AdamsMoultonIntegrator integrator = new AdamsMoultonIntegrator(2, h*1.0e-4, h, 1.0e-3, 1.0e-2);
        //GraggBulirschStoerIntegrator integrator = new GraggBulirschStoerIntegrator(1.0e-18, h, 1.0e-3, 1.0e-2);
        //DormandPrince54Integrator integrator = new DormandPrince54Integrator(1.0e-18, h, 1.0e-3, 1.0e-2);
        //ClassicalRungeKuttaIntegrator integrator = new ClassicalRungeKuttaIntegrator(h);
        integrator.addStepHandler(writeHandler);
        integrator.addStepHandler(term.getDetector());
        integrator.addEventHandler(term, Double.POSITIVE_INFINITY, 1.0e-12, 100);

        System.out.println("w: " + w);
        integrator.integrate(odes, 0, y0, tmax, y);
        
        System.out.println("final corr: " + SynchUtils.compCorr(y, n));
        
        if(correlate) {
            compCorr(params, y, dir + "two_time.txt");
        }

        long endTime = System.nanoTime();
        
        DynaComplexODEAdapter.toComplex(y, z0);
        PrintWriter writer = new PrintWriter(dir + "final_answer.txt", "UTF-8");
        for(int i = 0; i < dim; ++i) {
            writer.write(z0[i].getReal() + ", " + z0[i].getImaginary() + "\n");
        }
        writer.close();

        System.out.println("Run time: " + (endTime - startTime)/1.0e9 + " seconds");
    }

    private static double[] compCorr(CumulantParams params, double[] y, String filename) throws FileNotFoundException, UnsupportedEncodingException
    {
        int n = params.getN();
        
        int dim = SynchUtils.getDimension(n);
        int[] startIdx = new int[6];
        SynchUtils.getStartIdx(startIdx, n);
        DynaComplex[] z = new DynaComplex[dim];
        for(int i = 0; i < dim; ++i) {
            z[i] = new DynaComplex(0, 0);
        }
        DynaComplexODEAdapter.toComplex(y, z);
        
        DynaComplex[] szs = new DynaComplex[n];
        int idx2 = 0;
        for(int i = startIdx[2]; i < startIdx[2] + n; ++i) {
            szs[idx2] = new DynaComplex(z[i].getReal(), z[i].getImaginary());
            ++idx2;
        }
        
        double mod_sum = 0.0;
        DynaComplex[] z02 = new DynaComplex[n*n];
        idx2 = startIdx[3];
        for(int i = 0; i < n; ++i) {
            z02[i*n + i] = new DynaComplex(1.0, 0);
            for(int j = i+1; j < n; ++j) {
                z02[i*n + j] = new DynaComplex(z[idx2].getReal(), z[idx2].getImaginary());
                z02[j*n + i] = new DynaComplex(z[idx2].getReal(), z[idx2].getImaginary());
                z02[j*n + i].conjugate();
                
                mod_sum += 2.0*z02[j*n + i].getReal();
                
                //System.out.println(z02[i*n + j].getReal() + "+j*" + z02[i*n + j].getImaginary() + " " + idx2 + " " + startIdx[3] + " " + startIdx[4]);
                ++idx2;
            }
        }
        System.out.println("Avg corr factor: " + mod_sum/(n*(n-1)));

        DynaConstCoupling coupling = new DynaConstCoupling(params.getAlpha().getReal(), params.getAlpha().getImaginary());
        CorrelationODEs c_corr_odes = new CorrelationODEs(n, params.getGamma(), params.getW(), coupling, params.getD(), szs);
        DynaComplexODEAdapter odes = new DynaComplexODEAdapter(c_corr_odes);
        
        AdamsMoultonIntegrator integrator = new AdamsMoultonIntegrator(2, 1.0e-18, .001, 1.0e-3, 1.0e-2);
        if(!filename.isEmpty()) {
//            int[] out_col = {0,1,2,3,2*(n+1),2*(n+1)+1};
//            WriteHandler writeHandler = new WriteHandler(filename, out_col);
            TwoTimeHandler writeHandler = new TwoTimeHandler(filename);
            integrator.addStepHandler(writeHandler);
        }
        
        double[] y02 = new double[2*n*n];
        DynaComplexODEAdapter.toReal(z02, y02);
        
        double[] y2 = new double[2*n*n];
        
        double startTime = System.nanoTime();
        
        integrator.integrate(odes, 0, y02, 5, y2);
        
        double endTime = System.nanoTime();
        
        System.out.println("Correlation time: " + (endTime - startTime)/1.0e9 + " seconds");
        
        return y2;
    }
}
