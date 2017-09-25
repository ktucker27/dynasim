package ode;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

import utils.DynaComplex;

public class RPAAllToAllODEs implements FirstOrderDifferentialEquations {

    private int n;
    private double gamma;
    private double w;
    private double f;
    private double g;
    private double[] d;
    
    /**
     * The following convention is used for all Pauli operator types in this class:
     * 0 = x, 1 = y, 2 = z  
     */
    
    /**
     * vsums[i] = the sum of all <s_a^i>
     */
    double[] vsums;
    
    /**
     * rowsums[i][j][a] = sum_b <s_a^i s_b^j>
     */
    double[][][] rowsums;
    
    /**
     * Stores the start indices for each variable establishing the following order:
     * 0: <s_a^x>
     * 1: <s_a^y>
     * 2: <s_a^z>
     * 3: <s_a^x s_b^x>
     * 4: <s_a^y s_b^y>
     * 5: <s_a^z s_b^z>
     * 6: <s_a^x s_b^y>
     * 7: <s_a^x s_b^z>
     * 8: <s_a^y s_b^z>
     */
    int[] startIdx;
    
    /**
     * @param params input object specifying parameters
     */
    public RPAAllToAllODEs(CumulantParams params) {
        super();
        this.n = params.getN();
        this.gamma = params.getGamma();
        this.w = params.getW();
        this.f = params.getAlpha().getReal();
        this.g = params.getAlpha().getImaginary();
        this.d = params.getD();
        
        init();
    }
    
    /**
     * @param n number of atoms
     * @param gamma spontaneous decay rate
     * @param w pumping rate
     * @param feff effective mean field coupling
     * @param d size n vector of natural frequencies
     */
    public RPAAllToAllODEs(int n, double gamma, double w, double gel, DynaComplex alpha, double[] d) {
        super();
        this.n = n;
        this.gamma = gamma;
        this.w = w;
        this.f = alpha.getReal();
        this.g = alpha.getImaginary();
        this.d = d;
        
        init();
    }
    
    private void init() {
        startIdx = new int[9];
        startIdx[0] = 0;
        startIdx[1] = n;
        startIdx[2] = startIdx[1] + n;
        startIdx[3] = startIdx[2] + n;
        startIdx[4] = startIdx[3] + n*(n-1)/2;
        startIdx[5] = startIdx[4] + n*(n-1)/2;
        startIdx[6] = startIdx[5] + n*(n-1)/2;
        startIdx[7] = startIdx[6] + n*(n-1);
        startIdx[8] = startIdx[7] + n*(n-1);
        
        vsums = new double[3];
        rowsums = new double[3][3][n];
    }
    
    public void setW(double w) {
        this.w = w;
    }
    
    private int getTriIdx(int i, int j) {
        if(i > j) {
            int t = j;
            j = i;
            i = t;
        }
        
        return n*i - i*(i+1)/2 + j - i - 1;
    }
    
    private int getRecIdx(int i, int j) {
        if(j < i) {
            return i*(n-1) + j;
        }
        
        return i*(n-1) + j - 1;
    }
    
    /**
     * PRECONDITION: i < j
     * 
     * @param i First superscript
     * @param j Second superscript
     * @return Index offset for rectangular data array
     */
    private int getSuperTriIdx(int i, int j) {
        return 3*i - i*(i+1)/2 + j - i - 1;
    }
    
    private double getDouble(int al, int bt, int a, int b, double[] y) {
        if(al == bt) {
            return y[startIdx[3 + al] + getTriIdx(a,b)];
        }
        
        if(al < bt) {
            return y[startIdx[6 + getSuperTriIdx(al, bt)] + getRecIdx(a,b)];
        }
        
        return y[startIdx[6 + getSuperTriIdx(bt, al)] + getRecIdx(b,a)];
    }
    
    private double rpaSum(int al, int bt, int gm, int a, int b, int sidx, double[] y) {
        switch(sidx) {
        case 0:
            return y[startIdx[al] + a]*(rowsums[bt][gm][b] - getDouble(bt, gm, b, a, y));
        case 1:
            return y[startIdx[bt] + b]*(rowsums[al][gm][a] - getDouble(al, gm, a, b, y));
        case 2:
            return (vsums[gm] - y[startIdx[gm] + a] - y[startIdx[gm] + b])*getDouble(al, bt, a, b, y);
        }
        
        // This should never happen
        System.err.println("Invalid sidx passed to rpaSum: " + sidx);
        
        return 0.0;
    }
    
    @Override
    public void computeDerivatives(double t, double[] y, double[] yDot)
            throws MaxCountExceededException, DimensionMismatchException {
        
        // Start by computing the sums that we need.  This is O(n^2)

        // Zero out row and column sums
        for(int i = 0; i < 3; ++i) {
            vsums[i] = 0;
            for(int j = 0; j < 3; ++j) {
                for(int k = 0; k < n; ++k) {
                    rowsums[i][j][k] = 0;
                }
            }
        }

        for(int a = 0; a < n; ++a) {
            for(int i = 0; i < 3; ++i) {
                vsums[i] += y[startIdx[i] + a];
            }

            for(int b = 0; b < n; ++b) {
                if(a == b) continue;
                
                for(int i = 0; i < 3; ++i) {
                    for(int j = 0; j < 3; ++j) {
                        rowsums[i][j][a] += getDouble(i, j, a, b, y);
                    }
                }
            }
        }
        
        for(int a = 0; a < n; ++a) {
            // s^x
            yDot[a] = 0.5*gamma*(f*(rowsums[2][2][a] + rowsums[1][1][a]) - g*(rowsums[1][2][a] - rowsums[2][1][a])) 
                    - (gamma + w)*y[a] - (w - gamma);
        
            // s^y
            yDot[startIdx[1] + a] = -1.0*d[a]*y[startIdx[2] + a] - 0.5*(gamma + w)*y[startIdx[1] + a]
                    - 0.5*gamma*(f*rowsums[0][1][a] - g*rowsums[0][2][a]);
            
            // s^z
            yDot[startIdx[2] + a] = d[a]*y[startIdx[1] + a] - 0.5*(gamma + w)*y[startIdx[2] + a]
                    - 0.5*gamma*(f*rowsums[0][2][a] + g*rowsums[0][1][a]);
        }
        
        int idx = 0;
        for(int a = 0; a < n; ++a) {
            for(int b = a + 1; b < n; ++b) {
                // s^x s^x
                yDot[startIdx[3] + idx] = 0.5*gamma*(f*(rpaSum(2, 0, 2, a, b, 2, y) - rpaSum(1, 0, 1, a, b, 0, y) + rpaSum(0, 2, 2, a, b, 2, y) - rpaSum(0, 1, 1, a, b, 0, y)) 
                                                   - g*(rpaSum(1, 0, 2, a, b, 2, y) - rpaSum(2, 0, 1, a, b, 0, y) + rpaSum(0, 1, 2, a, b, 2, y) - rpaSum(0, 2, 1, a, b, 1, y)))
                                          + (gamma - w)*(y[a] + y[b]) - 2*(gamma + w)*getDouble(0, 0, a, b, y)
                                          + gamma*f*(getDouble(2, 2, a, b, y) + getDouble(1, 1, a, b, y));
                
                // s^y s^y
                yDot[startIdx[4] + idx] = -1.0*d[b]*getDouble(1, 2, a, b, y) + gamma*f*(getDouble(0, 0, a, b, y) - 0.5*(y[a] + y[b]))
                                          - (gamma + w)*getDouble(1, 1, a, b, y) - d[a]*getDouble(2, 1, a, b, y)
                                          + 0.5*gamma*(f*(rpaSum(0, 1, 1, a, b, 0, y) + rpaSum(1, 0, 1, a, b, 0, y))
                                                     + g*(rpaSum(0, 1, 2, a, b, 2, y) + rpaSum(1, 0, 2, a, b, 2, y)));
                
                // s^z s^z
                yDot[startIdx[5] + idx] = d[a]*getDouble(1, 2, a, b, y) + gamma*f*(getDouble(0, 0, a, b, y) - 0.5*(y[a] + y[b]))
                        - (gamma + w)*getDouble(2, 2, a, b, y) + d[b]*getDouble(2, 1, a, b, y)
                        + 0.5*gamma*(f*(rpaSum(0, 2, 2, a, b, 2, y) + rpaSum(2, 0, 2, a, b, 2, y))
                                   + g*(rpaSum(0, 2, 1, a, b, 1, y) + rpaSum(2, 0, 1, a, b, 0, y)));
                
                ++idx;
            }
        }
        
        idx = 0;
        for(int a = 0; a < n; ++a) {
            for(int b = 0; b < n; ++b) {
                if(a == b) continue;
                
                // s^x s^y
                yDot[startIdx[6] + idx] = -1.0*d[b]*getDouble(0, 2, a, b, y) - 1.5*(gamma + w)*getDouble(0, 1, a, b, y) + (gamma - w)*y[startIdx[1] + b]
                                          + 0.5*gamma*(f*y[startIdx[1] + a] + g*y[startIdx[2] + a]) - gamma*f*getDouble(1, 0, a, b, y)
                                          - 0.5*gamma*(f*rpaSum(0, 0, 1, a, b, 0, y) - g*rpaSum(0, 0, 2, a, b, 2, y))
                                          + 0.5*gamma*(f*(rpaSum(2, 1, 2, a, b, 2, y) + rpaSum(1, 1, 1, a, b, 0, y))
                                                     - g*(rpaSum(1, 1, 2, a, b, 2, y) - rpaSum(2, 1, 1, a, b, 0, y)));

                // s^x s^z
                yDot[startIdx[7] + idx] = d[b]*getDouble(0, 1, a, b, y) - 1.5*(gamma + w)*getDouble(0, 2, a, b, y) + (gamma - w)*y[startIdx[2] + b]
                                          + 0.5*gamma*(f*y[startIdx[2] + a] - g*y[startIdx[1] + a]) - gamma*f*getDouble(2, 0, a, b, y)
                                          + 0.5*gamma*(f*rpaSum(0, 0, 2, a, b, 2, y) - g*rpaSum(0, 0, 1, a, b, 0, y))
                                          + 0.5*gamma*(f*(rpaSum(2, 2, 2, a, b, 2, y) + rpaSum(1, 2, 1, a, b, 1, y))
                                                     - g*(rpaSum(1, 2, 2, a, b, 2, y) - rpaSum(2, 2, 1, a, b, 0, y)));

                // s^y s^z
                yDot[startIdx[8] + idx] = -1.0*d[a]*getDouble(2, 2, a, b, y) + 0.5*gamma*g*(y[a] - y[b]) 
                                         - (gamma + w)*getDouble(1, 2, a, b, y) + d[b]*getDouble(1, 1, a, b, y)
                                         - 0.5*gamma*(f*(rpaSum(0, 2, 1, a, b, 1, y) + rpaSum(1, 0, 2, a, b, 2, y))
                                                    + g*(rpaSum(1, 0, 1, a, b, 0, y) - rpaSum(0, 2, 2, a, b, 2, y)));
                
                ++idx;
            }
        }
    }

    @Override
    public int getDimension() {
        return 3*n + 3*n*(n-1)/2 + 3*n*(n-1);
    }
    
    public int getStartIdx(int idx) {
        return startIdx[idx];
    }
}
