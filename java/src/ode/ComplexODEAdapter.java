package ode;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

public class ComplexODEAdapter implements FirstOrderDifferentialEquations {

    ComplexODEs myOdes;
    Complex[] z;
    Complex[] zDot;
    int n;
    
    public ComplexODEAdapter(ComplexODEs odes) {
        myOdes = odes;
        
        n = odes.getDimension();
        z = new Complex[n];
        zDot = new Complex[n];
    }
    
    @Override
    public int getDimension() {
        return myOdes.getDimension()*2;
    }

    @Override
    public void computeDerivatives(double t, double[] y, double[] yDot)
            throws MaxCountExceededException, DimensionMismatchException {
        // Form complex dependent variables
        toComplex(y, z);
        
        // Calculate complex derivative values
        myOdes.computeDerivatives(t, z, zDot);
        
        // Extract real and imaginary parts
        toReal(zDot, yDot);
    }

    public static void toComplex(double[] y, Complex[] z) {
        int m = y.length/2;
        for(int i = 0; i < m; ++i) {
            z[i] = new Complex(y[2*i], y[2*i+1]);
        }
    }
    
    public static void toReal(Complex[] z, double[] y) {
        int m = z.length;
        for(int i = 0; i < m; ++i) {
            y[2*i] = z[i].getReal();
            y[2*i+1] = z[i].getImaginary();
        }
    }
}
