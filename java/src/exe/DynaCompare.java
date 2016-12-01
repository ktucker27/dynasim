package exe;

import java.io.FileNotFoundException;
import java.io.UnsupportedEncodingException;

import ode.CorrelationODEs;
import ode.CumulantAllToAllODEs;
import ode.DynaComplexODEAdapter;

import org.apache.commons.math3.distribution.TDistribution;
import org.apache.commons.math3.ode.nonstiff.AdamsMoultonIntegrator;
import org.apache.commons.math3.random.Well19937c;

import utils.DynaComplex;
import utils.WriteHandler;
import coupling.DynaConstCoupling;

public class DynaCompare {

    /**
     * @param args
     * @throws FileNotFoundException 
     * @throws UnsupportedEncodingException 
     */
    public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException {
        
        int n = 512;
        double h = 0.01;
        double w = 2.5;
        double gamma = 1.0;
        boolean correlate = false;
        
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
//        double delta = 2.5;
//        TDistribution tdist = new TDistribution(new Well19937c(1), 1);
//        double[] d = tdist.sample(n);
//        for(int i = 0; i < n; ++i) {
//            d[i] = d[i]*delta;
//            //System.out.println(d[i]);
//        }
        
        // Set natural frequencies to zero
        double[] d = new double[n];
        for(int i = 0; i < n; ++i) {
            d[i] = 0.0;
        }
        
        DynaComplex alpha = new DynaComplex(0.2, 0.0);
        DynaConstCoupling coupling = new DynaConstCoupling(1.0, 2.0);
        //LinearCoupling coupling = new LinearCoupling(2.0, 0.5, n);
        
        //DynaCumulantODEs codes = new DynaCumulantODEs(n, gamma, w, coupling, d);
        CumulantAllToAllODEs codes = new CumulantAllToAllODEs(n, gamma, w, alpha, d);
        DynaComplexODEAdapter odes = new DynaComplexODEAdapter(codes);
        int dim = codes.getDimension();
        
        DynaComplex[] z0 = new DynaComplex[dim];
        for(int i = 0; i < dim; ++i) {
            z0[i] = new DynaComplex(0, 0);
        }
        
        for(int i = 0; i < n; ++i) {
            z0[i] = new DynaComplex(0.2, 0.0);
        }
        
        double[] y0 = new double[2*dim];
        DynaComplexODEAdapter.toReal(z0, y0);
        
        int[] out_col = {0, 1, 2, 3,
                        codes.getStartIdx(1), codes.getStartIdx(1) + 1, codes.getStartIdx(1) + 2, codes.getStartIdx(1) + 3, 
                        codes.getStartIdx(2), codes.getStartIdx(2) + 1, codes.getStartIdx(2) + 2, codes.getStartIdx(2) + 3,
                        codes.getStartIdx(3), codes.getStartIdx(3) + 1,
                        codes.getStartIdx(4), codes.getStartIdx(4) + 1,
                        codes.getStartIdx(5), codes.getStartIdx(5) + 1};
        
        WriteHandler writeHandler = new WriteHandler("/Users/kristophertucker/output/sim_N512_D0_g0_ic0p2_2.txt", out_col);
        AdamsMoultonIntegrator integrator = new AdamsMoultonIntegrator(2, h*1.0e-2, h, 1.0e-3, 1.0e-2);
        //GraggBulirschStoerIntegrator integrator = new GraggBulirschStoerIntegrator(1.0e-18, h, 1.0e-3, 1.0e-2);
        //DormandPrince54Integrator integrator = new DormandPrince54Integrator(1.0e-18, h, 1.0e-3, 1.0e-2);
        //ClassicalRungeKuttaIntegrator integrator = new ClassicalRungeKuttaIntegrator(h);
        integrator.addStepHandler(writeHandler);
        
        double[] y = new double[2*dim];

        long startTime = System.nanoTime();

        integrator.integrate(odes, 0, y0, 10, y);
        
        long endTime = System.nanoTime();
        
        System.out.println("Run time: " + (endTime - startTime)/1.0e9 + " seconds");
        
        // Compute the correlation function if requested
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
            
            writeHandler = new WriteHandler("/Users/kristophertucker/output/correlation/corr_N32_al2p0_D7p5.txt", out_col2);
            integrator = new AdamsMoultonIntegrator(2, 1.0e-18, .001, 1.0e-3, 1.0e-2);
            integrator.clearStepHandlers();
            integrator.addStepHandler(writeHandler);
            
            double[] y02 = new double[2*n*n];
            DynaComplexODEAdapter.toReal(z02, y02);
            
            double[] y2 = new double[2*n*n];
            
            startTime = System.nanoTime();
            
            integrator.integrate(odes, 0, y02, 5, y2);
            
            endTime = System.nanoTime();
            
            System.out.println("Correlation time: " + (endTime - startTime)/1.0e9 + " seconds");
        }
    }

}
