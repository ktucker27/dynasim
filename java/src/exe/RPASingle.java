package exe;

import java.io.FileNotFoundException;
import java.io.UnsupportedEncodingException;

import org.apache.commons.math3.ode.nonstiff.AdamsMoultonIntegrator;

import handlers.WriteHandlerRPA;
import ode.CumulantParams;
import ode.RPAAllToAllODEs;
import utils.DynaComplex;
import utils.SynchUtils;

public class RPASingle {

    /**
     * @param args
     * @throws UnsupportedEncodingException 
     * @throws FileNotFoundException 
     */
    public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException {
        int n = 4;
        double h = 0.001;
        double gamma = 1.0;
        double tmax = 5.0;
        //double tmin = 5.0;
        double delta = 5.0;
        double f = 1;
        double g = 5.0;
        
        double w = SynchUtils.getWOpt(n);
        w = 17.0;
        
        // Cooper
        boolean cooper = false;
        if(cooper) {
            n = 2;
            delta = 0;
            f = 0;
            w = 0;
            g = 4*9*2*Math.PI*52;
//            g = 50;
        }
        
        System.out.println("w = " + w);
        
        DynaComplex alpha = new DynaComplex(f, g);
        
        double[] d = new double[n];
        SynchUtils.detuneGauss(delta, d);
        
        CumulantParams params = new CumulantParams(n, gamma, w, delta, alpha, d);
        
        RPAAllToAllODEs odes = new RPAAllToAllODEs(params);
        
        // Initialize everyone to spin-up along the z-direction
        int dim = odes.getDimension();
        double[] y0 = new double[dim];
        for(int i = 0; i < dim; ++i) {
            y0[i] = 0.0;
        }
        
        int idx = 0;
        for(int i = 0; i < n; ++i) {
            y0[odes.getStartIdx(2) + i] = 1.0;
            for(int j = i + 1; j < n; ++j) {
                y0[odes.getStartIdx(5) + idx] = 1.0;
                ++idx;
            }
        }
        
        System.out.println(params.toString());
        
        int[] out_col = {odes.getStartIdx(0), odes.getStartIdx(0) + 1,
                         odes.getStartIdx(1), odes.getStartIdx(1) + 1,
                         odes.getStartIdx(2), odes.getStartIdx(2) + 1,
                         odes.getStartIdx(3), odes.getStartIdx(4), odes.getStartIdx(5),
                         odes.getStartIdx(6), odes.getStartIdx(7), odes.getStartIdx(8)};
        
//        WriteHandler writeHandler = new WriteHandler("/Users/tuckerkj/output/temp/rpa_" + params.getFilename(), out_col, false);
        WriteHandlerRPA writeHandler = new WriteHandlerRPA("/Users/kristophertucker/output/temp/rpa_corr_" + params.getFilename(), params.getN());
        AdamsMoultonIntegrator integrator = new AdamsMoultonIntegrator(2, h*1.0e-4, h, 1.0e-3, 1.0e-2);
        integrator.addStepHandler(writeHandler);
        
        double[] y = new double[dim];
        
        long startTime = System.nanoTime();
        integrator.integrate(odes, 0, y0, tmax, y);
        
        long endTime = System.nanoTime();
        System.out.println("Run time: " + (endTime - startTime)/1.0e9 + " seconds");
    }

}
