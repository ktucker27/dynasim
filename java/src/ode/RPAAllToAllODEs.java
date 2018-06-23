package ode;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

import utils.DynaComplex;
import eval.RPAEval;

public class RPAAllToAllODEs implements FirstOrderDifferentialEquations {

    private int n;
    private double gamma;
    private double w;
    private double f;
    private double g;
    private double[] d;
    
    RPAEval myEval;
    
    /**
     * @param params input object specifying parameters
     */
    public RPAAllToAllODEs(SystemParams params) {
        super();
        this.n = params.getN();
        this.gamma = params.getGamma();
        this.w = params.getW();
        this.f = params.getAlpha().getReal();
        this.g = params.getAlpha().getImaginary();
        this.d = params.getD();
        
        myEval = new RPAEval(this.n);
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
        
        myEval = new RPAEval(this.n);
    }
    
    public void setW(double w) {
        this.w = w;
    }
    
    
    
    @Override
    public void computeDerivatives(double t, double[] y, double[] yDot)
            throws MaxCountExceededException, DimensionMismatchException {
        
        myEval.setVals(y);
        
        double[][][] rowsums = myEval.getRowSums();
        int[] startIdx = myEval.getStartIdx();
        
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
                yDot[startIdx[3] + idx] = 0.5*gamma*(f*(myEval.evalTriple(2, 0, 2, a, b) + myEval.evalTriple(1, 0, 1, a, b) + myEval.evalTriple(0, 2, 2, a, b) + myEval.evalTriple(0, 1, 1, a, b)) 
                                                   - g*(myEval.evalTriple(1, 0, 2, a, b) - myEval.evalTriple(2, 0, 1, a, b) + myEval.evalTriple(0, 1, 2, a, b) - myEval.evalTriple(0, 2, 1, a, b)))
                                          + (gamma - w)*(y[a] + y[b]) - 2*(gamma + w)*myEval.getDouble(0, 0, a, b, y)
                                          + gamma*f*(myEval.getDouble(2, 2, a, b, y) + myEval.getDouble(1, 1, a, b, y));
                
                // s^y s^y
                yDot[startIdx[4] + idx] = -1.0*d[b]*myEval.getDouble(1, 2, a, b, y) + gamma*f*(myEval.getDouble(0, 0, a, b, y) - 0.5*(y[a] + y[b]))
                                          - (gamma + w)*myEval.getDouble(1, 1, a, b, y) - d[a]*myEval.getDouble(2, 1, a, b, y)
                                          - 0.5*gamma*(f*(myEval.evalTriple(0, 1, 1, a, b) + myEval.evalTriple(1, 0, 1, a, b))
                                                     - g*(myEval.evalTriple(0, 1, 2, a, b) + myEval.evalTriple(1, 0, 2, a, b)));
                
                // s^z s^z
                yDot[startIdx[5] + idx] = d[a]*myEval.getDouble(1, 2, a, b, y) + gamma*f*(myEval.getDouble(0, 0, a, b, y) - 0.5*(y[a] + y[b]))
                        - (gamma + w)*myEval.getDouble(2, 2, a, b, y) + d[b]*myEval.getDouble(2, 1, a, b, y)
                        - 0.5*gamma*(f*(myEval.evalTriple(0, 2, 2, a, b) + myEval.evalTriple(2, 0, 2, a, b))
                                   + g*(myEval.evalTriple(0, 2, 1, a, b) + myEval.evalTriple(2, 0, 1, a, b)));
                
                ++idx;
            }
        }
        
        idx = 0;
        for(int a = 0; a < n; ++a) {
            for(int b = 0; b < n; ++b) {
                if(a == b) continue;
                
                // s^x s^y
                yDot[startIdx[6] + idx] = -1.0*d[b]*myEval.getDouble(0, 2, a, b, y) - 1.5*(gamma + w)*myEval.getDouble(0, 1, a, b, y) + (gamma - w)*y[startIdx[1] + b]
                                          + 0.5*gamma*(f*y[startIdx[1] + a] + g*y[startIdx[2] + a]) - gamma*f*myEval.getDouble(1, 0, a, b, y)
                                          - 0.5*gamma*(f*myEval.evalTriple(0, 0, 1, a, b) - g*myEval.evalTriple(0, 0, 2, a, b))
                                          + 0.5*gamma*(f*(myEval.evalTriple(2, 1, 2, a, b) + myEval.evalTriple(1, 1, 1, a, b))
                                                     - g*(myEval.evalTriple(1, 1, 2, a, b) - myEval.evalTriple(2, 1, 1, a, b)));

                // s^x s^z
                yDot[startIdx[7] + idx] = d[b]*myEval.getDouble(0, 1, a, b, y) - 1.5*(gamma + w)*myEval.getDouble(0, 2, a, b, y) + (gamma - w)*y[startIdx[2] + b]
                                          + 0.5*gamma*(f*y[startIdx[2] + a] - g*y[startIdx[1] + a]) - gamma*f*myEval.getDouble(2, 0, a, b, y)
                                          - 0.5*gamma*(f*myEval.evalTriple(0, 0, 2, a, b) + g*myEval.evalTriple(0, 0, 1, a, b))
                                          + 0.5*gamma*(f*(myEval.evalTriple(2, 2, 2, a, b) + myEval.evalTriple(1, 2, 1, a, b))
                                                     - g*(myEval.evalTriple(1, 2, 2, a, b) - myEval.evalTriple(2, 2, 1, a, b)));

                // s^y s^z
                yDot[startIdx[8] + idx] = -1.0*d[a]*myEval.getDouble(2, 2, a, b, y) + 0.5*gamma*g*(y[a] - y[b]) 
                                         - (gamma + w)*myEval.getDouble(1, 2, a, b, y) + d[b]*myEval.getDouble(1, 1, a, b, y)
                                         - 0.5*gamma*(f*(myEval.evalTriple(0, 2, 1, a, b) + myEval.evalTriple(1, 0, 2, a, b))
                                                    + g*(myEval.evalTriple(1, 0, 1, a, b) - myEval.evalTriple(0, 2, 2, a, b)));
                
                ++idx;
            }
        }
    }

    @Override
    public int getDimension() {
        return myEval.getDimension();
    }
    
    public int getStartIdx(int idx) {
        return myEval.getStartIdx()[idx];
    }
}
