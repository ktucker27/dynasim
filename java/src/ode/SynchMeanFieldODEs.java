package ode;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import java.lang.Math;

public class SynchMeanFieldODEs implements FirstOrderDifferentialEquations {

    private int n;
    private double gamma;
    private double w;
    private double feff;
    private double geff;
    private double[] d;

    /**
     * @param n number of oscillators
     * @param gamma spontaneous decay rate
     * @param w pumping rate
     * @param feff effective mean field coupling
     * @param d size n vector of natural frequencies
     */
    public SynchMeanFieldODEs(int n, double gamma, double w, double feff, double geff, double[] d) {
        super();
        this.n = n;
        this.gamma = gamma;
        this.w = w;
        this.feff = feff;
        this.geff = geff;
        this.d = d;
    }

    /**
     * @param t the independent time variable (not used)
     * @param y the 3*n dependent variables in the order (sz, sp, phi)
     * @param ydot the time derivative of y
     */
    @Override
    public void computeDerivatives(double t, double[] y, double[] ydot)
            throws MaxCountExceededException, DimensionMismatchException {
        // Check dimensions
        if(y.length != getDimension()) {
            throw new DimensionMismatchException(y.length, getDimension());
        }
        
        if(ydot.length != getDimension()) {
            throw new DimensionMismatchException(ydot.length, getDimension());
        }
        
        // Compute the order parameter
        Complex z = new Complex(0.0);
        for(int i = 0; i < n; ++i) {
            z = z.add(Complex.I.multiply(y[2*n+i]).exp().multiply(y[n+i]));
        }
        z = z.multiply(1/(double)n);
        
        double r = z.abs();
        double psi = z.getArgument();
        
        // Compute ydot
        for(int i = 0; i < n; ++i) {
            ydot[i] = r*y[n+i]*(-feff*Math.cos(psi - y[2*n+i]) + geff*Math.sin(psi - y[2*n+i])) - gamma*(0.5 + y[i]) + w*(0.5 - y[i]);
            ydot[n+i] = r*y[i]*(feff*Math.cos(psi - y[2*n+i]) - geff*Math.sin(psi - y[2*n+i])) - 0.5*(gamma + w)*y[n+i];
            ydot[2*n+i] = d[i] + r*((y[i]/y[n+i])*(feff*Math.sin(psi - y[2*n+i]) + geff*Math.cos(psi - y[2*n+i])));
        }
    }

    @Override
    public int getDimension() {
        return 3*n;
    }

}
