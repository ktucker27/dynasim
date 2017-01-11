package exe;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;

import ode.CumulantAllToAllODEs;
import ode.DynaComplexODEAdapter;

import org.apache.commons.math3.distribution.TDistribution;
import org.apache.commons.math3.ode.nonstiff.AdamsMoultonIntegrator;
import org.apache.commons.math3.random.Well19937c;

import utils.DynaComplex;
import utils.SynchUtils;

public class CumulantSearch {

    /**
     * @param args
     * @throws UnsupportedEncodingException 
     * @throws FileNotFoundException 
     */
    public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException {
        int n = 14;
        double h = 0.01;
        double wmin = 0.5;
        double wmax = 10.0;
        double dw = 0.5;
        double gamma = 1.0;
        
        // Get natural frequencies from Cauchy distribution
        double delta = 2.0;
        TDistribution tdist = new TDistribution(new Well19937c(1), 1);
        double[] d = tdist.sample(n);
        for(int i = 0; i < n; ++i) {
            d[i] = d[i]*delta;
            //System.out.println(d[i]);
        }
        
        // Set natural frequencies to zero
//        double[] d = new double[n];
//        for(int i = 0; i < n; ++i) {
//            d[i] = 1.0;
//        }
        
        DynaComplex alpha = new DynaComplex(1, 0.5);
        //DynaConstCoupling coupling = new DynaConstCoupling(1.0, 2.0);
        //LinearCoupling coupling = new LinearCoupling(2.0, 0.5, n);
        
        //DynaCumulantODEs codes = new DynaCumulantODEs(n, gamma, w, coupling, d);
        CumulantAllToAllODEs codes = new CumulantAllToAllODEs(n, gamma, 0, alpha, d);
        DynaComplexODEAdapter odes = new DynaComplexODEAdapter(codes);
        int dim = codes.getDimension();
        
        DynaComplex[] z0 = new DynaComplex[dim];
//        for(int i = 0; i < dim; ++i) {
//            z0[i] = new DynaComplex(0, 0);
//        }
//        
//        for(int i = 0; i < n; ++i) {
//            z0[i] = new DynaComplex(0.2, 0.0);
//        }
        
        SynchUtils.initialize(z0, Math.PI/2.0, n);
        
        double[] y0 = new double[2*dim];
        DynaComplexODEAdapter.toReal(z0, y0);
        
        AdamsMoultonIntegrator integrator = new AdamsMoultonIntegrator(2, h*1.0e-4, h, 1.0e-3, 1.0e-2);
        //GraggBulirschStoerIntegrator integrator = new GraggBulirschStoerIntegrator(1.0e-18, h, 1.0e-3, 1.0e-2);
        //DormandPrince54Integrator integrator = new DormandPrince54Integrator(1.0e-18, h, 1.0e-3, 1.0e-2);
        //ClassicalRungeKuttaIntegrator integrator = new ClassicalRungeKuttaIntegrator(h*1.0e-2);
        
        double[] y = new double[2*dim];

        long startTime = System.nanoTime();

        PrintWriter writer = new PrintWriter("/Users/kristophertucker/output/search_out.txt", "UTF-8");
        for(double w  = wmin; w <= wmax; w += dw) {
            System.out.println("w: " + w);
            codes.setW(w);
            integrator.integrate(odes, 0, y0, 10, y);

            DynaComplexODEAdapter.toComplex(y, z0);

            double sum = 0.0;
            for(int i = codes.getStartIdx(3)/2; i < codes.getStartIdx(3)/2 + n*(n-1)/2; ++i) {
                sum += 2.0*z0[i].getReal();
                //System.out.println(z0[i].toString());
            }
            sum *= 1.0/(n*(n-1));
            
            DynaComplex csum = new DynaComplex(0, 0);
            for(int i = 0; i < n; ++i) {
                csum.add(z0[i]);
            }
            csum.multiply(1.0/(double)n);
            
            writer.print(w + " " + csum.mod() + "\n");
        }
        writer.close();
        
        long endTime = System.nanoTime();

        System.out.println("Run time: " + (endTime - startTime)/1.0e9 + " seconds");
    }

}
