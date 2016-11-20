package ode;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

import utils.DynaComplex;

public class DynaComplexODEAdapter implements FirstOrderDifferentialEquations {

    DynaComplexODEs myOdes;
    DynaComplex[] z;
    DynaComplex[] zDot;
    int n;
    
    public DynaComplexODEAdapter(DynaComplexODEs odes) {
        myOdes = odes;
        
        n = odes.getDimension();
        z = new DynaComplex[n];
        zDot = new DynaComplex[n];
        for(int i = 0; i < n; ++i) {
            z[i] = new DynaComplex(0, 0);
            zDot[i] = new DynaComplex(0, 0);
        }
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
        
        // DEBUGG
//        System.out.println("------------------------------------------");
//        for(int i = 0; i < z.length; ++i) {
//            System.out.println(z[i].getReal() + " + j*" + z[i].getImaginary());
//        }
//        
//        System.out.println("");
//
//        for(int i = 0; i < z.length; ++i) {
//            System.out.println(zDot[i].getReal() + " + j*" + zDot[i].getImaginary());
//        }
        
        // Extract real and imaginary parts
        toReal(zDot, yDot);
    }

    public static void toComplex(double[] y, DynaComplex[] z) {
        int m = y.length/2;
        for(int i = 0; i < m; ++i) {
            z[i].setReal(y[2*i]);
            z[i].setImaginary(y[2*i+1]);
        }
    }
    
    public static void toReal(DynaComplex[] z, double[] y) {
        int m = z.length;
        for(int i = 0; i < m; ++i) {
            y[2*i] = z[i].getReal();
            y[2*i+1] = z[i].getImaginary();
        }
    }
}
