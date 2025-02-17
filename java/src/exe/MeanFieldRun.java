package exe;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.Random;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.ode.nonstiff.AdamsMoultonIntegrator;

import ode.SystemParams;
import ode.SynchMeanFieldODEs;
import utils.DynaComplex;
import utils.SynchUtils;

public class MeanFieldRun {

    /**
     * @param args
     * @throws UnsupportedEncodingException 
     * @throws FileNotFoundException 
     */
    public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException {
        long startTime = System.nanoTime();

        int n = 8;
        double delta = 1.0;
        double f = 1.0;
        double g = 30.0;
        
        double[] d = new double[n];
        SynchUtils.detuneGauss(delta, d);
//        SynchUtils.detuneLor(delta, d);
        
        double minw = 2.0;
        double maxw = 4.5;
        int numw = 25;
        double dw = (maxw - minw)/(double)(numw - 1);

        double[] r = new double[numw];
        double[] wv = new double[numw];
        double w = minw;
        for(int wi = 0; wi < numw; ++wi) {
            wv[wi] = w;
            r[wi] = runSim(w, f, g, d);
            System.out.println("w = " + w + " " + r[wi]);
            w += dw;
        }

        long endTime = System.nanoTime();

        SystemParams params = new SystemParams(n, 1, 0, delta, new DynaComplex(f, g), d);
//        PrintWriter writer = new PrintWriter("/Users/tuckerkj/output/mf/wrun_" + params.getFilename(), "UTF-8");
        PrintWriter writer = new PrintWriter("/Users/kristophertucker/output/mf/wrun_temp.txt", "UTF-8");
        for(int i = 0; i < numw; ++i) {
                writer.print(wv[i] + ", " + r[i] + "\n");
        }
        writer.close();

        System.out.println("Run time: " + (endTime - startTime)/1.0e9 + " seconds");
    }

    public static double runSim(double w, double feff, double geff, double[] d) {
        int n = d.length;
        double h = 0.001;

        SynchMeanFieldODEs odes = new SynchMeanFieldODEs(n, 1, w, feff, geff, d);
        
//        EulerIntegrator integrator = new EulerIntegrator(h);
        AdamsMoultonIntegrator integrator = new AdamsMoultonIntegrator(2, h*1.0e-4, h, 1.0e-3, 1.0e-2);

        double[] y = new double[3*n];
        double[] y0 = new double[3*n];
        
        Random rand = new Random(5);
        for(int i = 0; i < n; ++i) {
            y0[i] = 0.0;
            y0[n+i] = 1.0;
//            y0[2*n+i] = rand.nextDouble()*2*Math.PI - Math.PI;
            y0[2*n+i] = 0.0;
//            y0[2*n+i] = i*2*Math.PI/(double)n - Math.PI;
        }
        
        double t = 200.0;
        integrator.integrate(odes, 0, y0, t, y);
        
        // Compute the order parameter (|avg_a(<sigma_a^+>)|)
        Complex z = new Complex(0.0);
        for(int i = 0; i < n; ++i) {
            z = z.add(Complex.I.multiply(-y[2*n+i]).exp().multiply(y[n+i]));
        }
        z = z.multiply(1/(double)n);

        return z.abs();
        //double psi = z.getArgument();
    }
}
