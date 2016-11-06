package ode;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;

public interface ComplexODEs {
    public void computeDerivatives(double t, Complex[] z, Complex[] zDot)
            throws MaxCountExceededException, DimensionMismatchException;

    public int getDimension();
}
