package coupling;

import utils.DynaComplex;

/**
 * 
 * @author kristophertucker
 *
 * An interface for a complex coupling matrix.  The goal is to provide a matrix
 * like interface while allowing for implementations that do not have to store 
 * the entire matrix of values.  This allows us to exploit the symmetry that 
 * often comes with complex coupling matrices.
 */
public interface DynaComplexCoupling {
    double getF(int row, int col);

    double getG(int row, int col);

    DynaComplex getAlpha(int row, int col);
}
