package ode;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

import coupling.ComplexCoupling;

public class SynchODEs implements FirstOrderDifferentialEquations {
    private int n;
    private double gamma;
    private double w;
    private ComplexCoupling coupling;
    private double[] d;

    /**
     * @param n number of oscillators
     * @param gamma spontaneous decay rate
     * @param w pumping rate
     * @param feff effective mean field coupling
     * @param d size n vector of natural frequencies
     */
    public SynchODEs(int n, double gamma, double w, ComplexCoupling coupling, double[] d) {
        super();
        this.n = n;
        this.gamma = gamma;
        this.w = w;
        this.coupling = coupling;
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
        Complex[] z = new Complex[n];
        for(int k = 0; k < n; ++k) {
            z[k] = new Complex(0.0);
            for(int i = 0; i < n; ++i) {
                z[k] = z[k].add(Complex.I.multiply(y[2*n+i]).exp().multiply(y[n+i]).multiply(coupling.getAlpha(k, i)));
            }
        }
        
        // Compute ydot
        for(int i = 0; i < n; ++i) {
            Complex op = Complex.I.multiply(-y[2*n+i]).exp().multiply(z[i]);
            ydot[i] = -gamma*y[n+i]*op.getReal() - gamma*(0.5+y[i]) + w*(0.5-y[i]);
            ydot[n+i] = -0.5*(gamma + w)*y[n+i] + gamma*y[i]*op.getReal();
            ydot[2*n+i] = -d[i] + (gamma*y[i]/y[n+i])*op.getImaginary();
        }
    }

    @Override
    public int getDimension() {
        return 3*n;
    }
}
