package exe;

import java.io.FileNotFoundException;
import java.io.UnsupportedEncodingException;

import ode.ComplexODEAdapter;
import ode.ComplexTestODEs;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.ode.nonstiff.AdamsMoultonIntegrator;

import utils.WriteHandler;

public class ComplexTest {

    /**
     * @param args
     * @throws UnsupportedEncodingException 
     * @throws FileNotFoundException 
     */
    public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException {
        double w = 2*Math.PI;
        ComplexTestODEs odes = new ComplexTestODEs(w);
        
        ComplexODEAdapter adapter = new ComplexODEAdapter(odes);
        WriteHandler writeHandler = new WriteHandler("/Users/kristophertucker/Google Drive/Research/Synch/output/c_test.txt", new int[] {1});
        
        double h = .01;
        AdamsMoultonIntegrator integrator = new AdamsMoultonIntegrator(2, 1.0e-18, h, 1.0e-3, 1.0e-2);
        integrator.addStepHandler(writeHandler);
        
        double[] y = new double[2];
        double[] y0 = new double[2];
        Complex[] z = new Complex[1];
        z[0] = new Complex(1.0, 0.0);
        ComplexODEAdapter.toReal(z, y0);
        integrator.integrate(adapter, 0, y0, 10, y);
        
        ComplexODEAdapter.toComplex(y, z);
        
        System.out.println(z[0].getReal() + ", " + z[0].getImaginary());
        
        System.out.println("Done");
    }

}
