package exe;

import java.io.FileNotFoundException;
import java.io.UnsupportedEncodingException;
import java.util.Random;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.ode.nonstiff.AdamsMoultonIntegrator;

import handlers.RingSteadyStateTest;
import handlers.SynchSteadyStateTerminator;
import handlers.WriteHandlerMeanField;
import ode.SystemParams;
import ode.SynchMeanFieldODEs;
import utils.DynaComplex;
import utils.SynchUtils;

public class MeanFieldSingle {

    /**
     * @param args
     * @throws UnsupportedEncodingException 
     * @throws FileNotFoundException 
     */
    public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException {
        long startTime = System.nanoTime();

        int n = 8;
        double h = 0.001;
//        double a = 0.025;
        double delta = 1.0;
        double gamma = 1.0;
        double f = 1.0;
        double g = 30.0;
        double w = 3.042;
//      double tmin = 0.0;
        double tmax = 200.0;
        
        double[] d = new double[n];
        SynchUtils.detuneGauss(delta, d);
//        SynchUtils.detuneLor(delta, d);
        
        double[] y0 = new double[3*n];
        double[] y = new double[3*n];
        
        Random rand = new Random(5);
        for(int i = 0; i < n; ++i) {
            y0[i] = 0.0;
            y0[n+i] = 1.0;
//            y0[2*n+i] = rand.nextDouble()*2*Math.PI - Math.PI;
//            y0[2*n+i] = 0.0;
            y0[2*n+i] = i*2*Math.PI/(double)n - Math.PI;
        }
        
        SystemParams params = new SystemParams(n, gamma, w, delta, new DynaComplex(f,g), d);

        String dir = "/Users/tuckerkj/output/temp/";
//        WriteHandler writeHandler = new WriteHandler("/Users/kristophertucker/output/temp/.txt", new int[] {n+n/2}, true);
//        WriteAllHandler writeHandler = new WriteAllHandler("/Users/tuckerkj/output/temp/mf_all_N100_D8_g1_w10.txt", 3*n);
        WriteHandlerMeanField writeHandler = new WriteHandlerMeanField(dir + "mf_corr2_" + params.getFilename(), n);
        
        RingSteadyStateTest test = new RingSteadyStateTest();
        test.setQuietMode(false);
        SynchSteadyStateTerminator term = new SynchSteadyStateTerminator(test, 0.5, 500, 1500000);
        
        AdamsMoultonIntegrator integrator = new AdamsMoultonIntegrator(2, h*1.0e-4, h, 1.0e-3, 1.0e-2);
        integrator.addStepHandler(writeHandler);
//        integrator.addStepHandler(term.getDetector());
//        integrator.addEventHandler(term, Double.POSITIVE_INFINITY, 1.0e-12, 100);
        
        //RingCoupling coupling = new RingCoupling(n,a);
//        ConstCoupling coupling = new ConstCoupling(f, g);
//        SynchODEs odes = new SynchODEs(n, gamma, w, coupling, d);
        SynchMeanFieldODEs odes = new SynchMeanFieldODEs(n, gamma, w, f, g, d);
        try {
            integrator.integrate(odes, 0, y0, tmax, y);
//            System.out.println("R:");
//            for(int i = n; i < 2*n; ++i) {
//                System.out.print(y[i] + " ");
//            }
//            System.out.print("\n");

            // Compute the order parameter (|avg_a(<sigma_a^+>)|)
            Complex z = new Complex(0.0);
            for(int i = 0; i < n; ++i) {
                z = z.add(Complex.I.multiply(-y[2*n+i]).exp().multiply(y[n+i]));
            }
            z = z.multiply(1/(double)n);
            
            System.out.println("order: " + z.abs());
        } catch(Exception ex) {
            System.out.println(ex.getMessage());
        }

        long endTime = System.nanoTime();

        System.out.println("Run time: " + (endTime - startTime)/1.0e9 + " seconds");
    }
}
