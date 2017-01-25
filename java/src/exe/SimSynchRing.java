package exe;

import integrator.AdamsMoultonFactory;
import integrator.ThreadPoolIntegrator;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.Random;

import ode.SynchODEs;
import utils.FullODESolution;
import coupling.RingCoupling;

public class SimSynchRing {

    /**
     * @param args
     * @throws UnsupportedEncodingException 
     * @throws FileNotFoundException 
     */
    public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException {
        long startTime = System.nanoTime();

        int n = 103;
        double h = 0.0005;
        double a = 0.05;
        
        double[] d = new double[n];
        for(int i = 0; i < n; ++i) {
            d[i] = 0.0;
        }
        
        double[] y0 = new double[3*n];
        
        Random rand = new Random();
        for(int i = 0; i < n; ++i) {
            y0[i] = rand.nextDouble();
            y0[n+i] = rand.nextDouble();
            y0[2*n+i] = rand.nextDouble()*2*Math.PI - Math.PI;
        }

        double t0 = 0.0;
        double t = 10.0;
        
        double wStart = 1.0;
        double hw = 63.0/100.0;
        int numw = 25;
        
        //ThreadPoolIntegrator integrator = new ThreadPoolIntegrator(8, new EulerIntegratorFactory(h));
        ThreadPoolIntegrator integrator = new ThreadPoolIntegrator(8, new AdamsMoultonFactory(1.0e-18, h));
        integrator.setQuietMode(false);

        RingCoupling coupling = new RingCoupling(n,a);
        FullODESolution[] soln = new FullODESolution[numw];
        for(int j = 0; j < numw; ++j) {
            soln[j] = new FullODESolution();
            SynchODEs odes = new SynchODEs(n, 1, wStart + j*hw, coupling, d);
            integrator.addIvp(odes, t0, y0, t, soln[j]);
        }

        boolean success = integrator.waitForFinished(12*3600);

        long endTime = System.nanoTime();

        if(success) {
            PrintWriter writer = new PrintWriter("/Users/kristophertucker/Google Drive/Research/Synch/output/r_vs_w3.txt", "UTF-8");
            for(int j = 0; j < numw; ++j) {
                if(soln[j].getSolution() != null) {
                    writer.print(soln[j].getSolution()[n] + " ");
                    System.out.println(soln[j].getSolution()[n]);
                } else {
                    writer.print("0 ");
                }
            }
            writer.close();
        }

        System.out.println("Run time: " + (endTime - startTime)/1.0e9 + " seconds");
    }
}
