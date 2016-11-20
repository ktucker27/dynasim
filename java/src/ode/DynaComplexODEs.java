package ode;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;

import utils.DynaComplex;

public interface DynaComplexODEs {
    public void computeDerivatives(double t, DynaComplex[] z, DynaComplex[] zDot)
            throws MaxCountExceededException, DimensionMismatchException;

    public int getDimension();
}
