package utils;

import org.apache.commons.math3.complex.Complex;

public class OrderParameterSolution implements ODESolution {

    private double r;
    
    public OrderParameterSolution() {
        super();
        r = 0.0;
    }

    @Override
    public void setSolution(double[] y) {
        r = compOrderParam(y);
    }
    
    public double getOrderParam() {
        return r;
    }

    public static double compOrderParam(double[] y) {
        // Compute the order parameter
        int n = y.length/3;
        Complex z = new Complex(0.0);
        for(int i = 0; i < n; ++i) {
            z = z.add(Complex.I.multiply(y[2*n+i]).exp().multiply(y[n+i]));
        }
        z = z.multiply(1/(double)n);

        return z.abs();
        //double psi = z.getArgument();
    }
}
