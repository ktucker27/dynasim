package exe;

import java.io.FileNotFoundException;
import java.io.UnsupportedEncodingException;

import ode.DynaComplexODEAdapter;
import ode.DynaCumulantODEs;

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
        
        int n = 32;
        double h = 0.01;
        double delta = 7.5;
        double w = 2.5;
        double gamma = 1.0;
        
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
        TDistribution tdist = new TDistribution(new Well19937c(1), 1);
        double[] d = tdist.sample(n);
        for(int i = 0; i < n; ++i) {
            d[i] = d[i]*delta;
            System.out.println(d[i]);
        }
        
        DynaConstCoupling coupling = new DynaConstCoupling(1.0, 2.0);
        
        DynaCumulantODEs codes = new DynaCumulantODEs(n, gamma, w, coupling, d);
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
        
        WriteHandler writeHandler = new WriteHandler("/Users/kristophertucker/output/compare_hp01_N32_g2p0_ic0p2.txt", out_col);
        AdamsMoultonIntegrator integrator = new AdamsMoultonIntegrator(2, 1.0e-18, h, 1.0e-3, 1.0e-2);
        //GraggBulirschStoerIntegrator integrator = new GraggBulirschStoerIntegrator(1.0e-18, h, 1.0e-3, 1.0e-2);
        //DormandPrince54Integrator integrator = new DormandPrince54Integrator(1.0e-18, h, 1.0e-3, 1.0e-2);
        //ClassicalRungeKuttaIntegrator integrator = new ClassicalRungeKuttaIntegrator(h);
        integrator.addStepHandler(writeHandler);
        
        double[] y = new double[2*dim];

        long startTime = System.nanoTime();

        integrator.integrate(odes, 0, y0, 10, y);
        
        long endTime = System.nanoTime();
        
        System.out.println("Run time: " + (endTime - startTime)/1.0e9 + " seconds");
    }

}
