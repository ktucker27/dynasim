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
    
    public SynchMeanFieldODEs(SystemParams params) {
        super();
        this.n = params.getN();
        this.gamma = params.getGamma();
        this.w = params.getW();
        this.feff = params.getAlpha().getReal();
        this.geff = params.getAlpha().getImaginary();
        this.d = params.getD();
    }

    /**
     * @param t the independent time variable (not used)
     * @param y the 3*n dependent variables in the order (<sigma^z>, abs(<sigma^+>), arg(<sigma^+>))
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
            z = z.add(Complex.I.multiply(y[2*n+i]).exp().multiply(0.5*y[n+i]));
        }
        //z = z.multiply(1/(double)n);
        
        double r = z.abs();
        double psi = z.getArgument();
        
        // Compute ydot
        for(int i = 0; i < n; ++i) {
            ydot[i] = -0.5*y[n+i]*(feff*(2*r*Math.cos(psi - y[2*n+i]) - y[n+i]) + geff*2*r*Math.sin(psi - y[2*n+i])) - gamma*(1 + y[i]) + w*(1 - y[i]);
            ydot[n+i] = 0.5*gamma*y[i]*(feff*(2*r*Math.cos(psi - y[2*n+i]) - y[n+i]) + geff*2*r*Math.sin(psi - y[2*n+i])) - 0.5*(gamma + w)*y[n+i];
            ydot[2*n+i] = -d[i] + 0.5*gamma*((y[i]/y[n+i])*(feff*2*r*Math.sin(psi - y[2*n+i]) - geff*(2*r*Math.cos(psi - y[2*n+i]) - y[n+i])));
            
//            ydot[i] = -0.5*y[n+i]*(feff*(2*r*Math.cos(psi - y[2*n+i])) + geff*2*r*Math.sin(psi - y[2*n+i])) - gamma*(1 + y[i]) + w*(1 - y[i]);
//            ydot[n+i] = 0.5*gamma*y[i]*(feff*(2*r*Math.cos(psi - y[2*n+i])) + geff*2*r*Math.sin(psi - y[2*n+i])) - 0.5*(gamma + w)*y[n+i];
//            ydot[2*n+i] = -d[i] + 0.5*gamma*((y[i]/y[n+i])*(feff*2*r*Math.sin(psi - y[2*n+i]) - geff*(2*r*Math.cos(psi - y[2*n+i]))));
        }
    }

    @Override
    public int getDimension() {
        return 3*n;
    }

}
