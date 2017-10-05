package exe;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;

import ode.CumulantParams;
import ode.DynaComplexODEAdapter;
import ode.MasterAllToAllODEs;

import org.apache.commons.math3.ode.nonstiff.AdamsMoultonIntegrator;

import utils.DynaComplex;
import utils.SynchUtils;
import utils.TPSOperator;

public class MasterSingle {

    /**
     * @param args
     * @throws UnsupportedEncodingException 
     * @throws FileNotFoundException 
     */
    public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException {
        int n = 9;
        double h = 0.001;
        double gamma = 1.0;
        double tmax = 5.0;
        //double tmin = 5.0;
        double delta = 5.0;
        double f = 1;
        double g = 5.0;
        
        boolean outputZDist = true;
        
        double w = SynchUtils.getWOpt(n);
        w = 2.0;
        System.out.println("w = " + w);
        
        DynaComplex alpha = new DynaComplex(f, g);
        
        double[] d = new double[n];
        SynchUtils.detuneGauss(delta, d);
        
        CumulantParams params = new CumulantParams(n, gamma, w, delta, alpha, d);
        
        // Initialize everyone to spin-up along the x-direction
        TPSOperator rho0 = new TPSOperator(n);
        rho0.set(1.0/Math.pow(2,n));
        
        DynaComplex[] z0 = rho0.getVals();
        double[] y0 = new double[2*z0.length];
        DynaComplexODEAdapter.toReal(z0, y0);
        
        System.out.println(params.toString());
        
        MasterAllToAllODEs modes = new MasterAllToAllODEs(params);
        DynaComplexODEAdapter odes = new DynaComplexODEAdapter(modes);
        
        String dir = "/Users/kristophertucker/output/temp/";
        //WriteHandlerMaster writeHandler = new WriteHandlerMaster("/Users/kristophertucker/output/temp/master_" + params.getFilename(), n);
        AdamsMoultonIntegrator integrator = new AdamsMoultonIntegrator(2, h*1.0e-4, h, 1.0e-3, 1.0e-2);
        //integrator.addStepHandler(writeHandler);
        
        double[] y = new double[2*z0.length];
        
        long startTime = System.nanoTime();
        integrator.integrate(odes, 0, y0, tmax, y);
        
        long endTime = System.nanoTime();
        System.out.println("Run time: " + (endTime - startTime)/1.0e9 + " seconds");
        
        if(outputZDist) {
            DynaComplex[] z = new DynaComplex[z0.length];
            for(int i = 0; i < z.length; ++i) {
                z[i] = new DynaComplex();
            }
            DynaComplexODEAdapter.toComplex(y, z);
            TPSOperator rho = new TPSOperator(z);
            
            writeZDist(rho, dir + "zdist_" + params.getFilename());
        }
    }
    
    public static void writeZDist(TPSOperator rho, String filepath) throws FileNotFoundException {
        int n = rho.getN();
        
        double[] dist = new double[2*n + 1];
        for(int i = 0; i < dist.length; ++i) {
            dist[i] = 0.0;
        }
        
        for(int i = 0; i < SynchUtils.pow(2, n); ++i) {
            int k = 0;
            int ii = i;
            for(int j = 0; j < n; ++j) {
                k += (ii % 2);
                ii = ii >> 1;
            }
            int eval = n - 2*k;
            dist[eval + n] += rho.getVal(i, i).getReal();
        }
        
        PrintWriter writer = new PrintWriter(filepath);
        for(int i = 0; i < dist.length; ++i) {
            writer.write(dist[i] + "\n");
        }
        writer.close();
    }

}
