package exe;

import java.io.FileNotFoundException;
import java.io.UnsupportedEncodingException;

import mcwf.MCWFIntegrator;
import mcwf.MCWFWriter;
import ode.SystemParams;
import utils.DynaComplex;

public class MCWF {

    public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException {
        // Set simulation parameters
        int numTrajectories = 1;
        int numThreads = 1;
        
        int n = 4;
        double gaa = 2.0;
        double gab = 2.0;
        double o = 10.0;
        double w = 0.0;
        double faa = 0.0;
        double fab = 0.0;
        double gel = 0.0;
        double gamma = 1.0;
        
        double dt = 1.0e-4;
        double tf = 1.0;
        int numTimes = (int)Math.ceil(tf/dt + 1);
        
        long timeout = 7*24*3600; // One week
        
        double[] d = new double[n];
        for(int i = 0; i < n; ++i) {
            d[i] = 0.0;
        }
        
        SystemParams params = new SystemParams(n, gamma, w, o, gel, 0, faa, fab, gaa, gab, d);
        
        // Define the initial condition
        DynaComplex[] initialState = new DynaComplex[n+1];
        for(int i = 0; i < n+1; ++i) {
            initialState[i] = new DynaComplex(0,0);
        }
        initialState[n].set(1,0);
        
        // Create the integrator
        MCWFIntegrator integrator = new MCWFIntegrator(numTrajectories, numTimes, dt, params, initialState, numThreads);
        
        // Integrate
        long startTime = System.nanoTime();
        integrator.start();
        integrator.waitForFinished(timeout);

        // Write results
        MCWFWriter writer = new MCWFWriter("/Users/tuckerkj/output/mcwf_test.csv");
        writer.write(integrator.getAggregator());
        
        long endTime = System.nanoTime();
        System.out.println("Run time: " + (endTime - startTime)/1.0e9 + " seconds");
    }

}
