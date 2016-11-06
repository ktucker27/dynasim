package ode;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;

public class ComplexTestODEs implements ComplexODEs {

    double w;
    
    public ComplexTestODEs(double ww) {
        w = ww;
    }
    
    @Override
    public void computeDerivatives(double t, Complex[] y, Complex[] yDot)
            throws MaxCountExceededException, DimensionMismatchException {
        yDot[0] = new Complex(-w*y[0].getImaginary(), w*y[0].getReal());
    }

    @Override
    public int getDimension() {
        return 1;
    }

}
