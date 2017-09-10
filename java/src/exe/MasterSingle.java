package exe;

import java.io.FileNotFoundException;
import java.io.UnsupportedEncodingException;

import handlers.WriteHandlerMaster;

import org.apache.commons.math3.ode.nonstiff.AdamsMoultonIntegrator;

import ode.CumulantParams;
import ode.DynaComplexODEAdapter;
import ode.MasterAllToAllODEs;
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
        int n = 2;
        double h = 0.001;
        double gamma = 1.0;
        double tmax = 5.0;
        double tmin = 5.0;
        double delta = 0.0;
        double f = 1;
        double g = 5.0;
        
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
        
        WriteHandlerMaster writeHandler = new WriteHandlerMaster("/Users/kristophertucker/output/temp/master_out.txt", n);
        AdamsMoultonIntegrator integrator = new AdamsMoultonIntegrator(2, h*1.0e-4, h, 1.0e-3, 1.0e-2);
        integrator.addStepHandler(writeHandler);
        
        double[] y = new double[2*z0.length];
        
        long startTime = System.nanoTime();
        integrator.integrate(odes, 0, y0, tmax, y);
        
        long endTime = System.nanoTime();
        System.out.println("Run time: " + (endTime - startTime)/1.0e9 + " seconds");
    }

}
