package ode;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;

import utils.DynaComplex;
import utils.SynchUtils;
import utils.TPSOperator;
import utils.TPSOperator.PauliOp;

public class MasterAllToAllODEs implements DynaComplexODEs {

    private int n;
    private double gamma;
    private double w;
    private DynaComplex alpha;
    private double[] d;
    
    private DynaComplex t1;
    private DynaComplex t2;
    
    TPSOperator sum1;
    
    /**
     * @param params input object specifying parameters
     */
    public MasterAllToAllODEs(SystemParams params) {
        super();
        this.n = params.getN();
        this.gamma = params.getGamma();
        this.w = params.getW();
        this.alpha = new DynaComplex(params.getAlpha());
        this.d = params.getD();
        
        init();
    }
    
    private void init() {
        t1 = new DynaComplex();
        t2 = new DynaComplex();
        
        sum1 = new TPSOperator(n);
    }
    
    @Override
    public void computeDerivatives(double t, DynaComplex[] z, DynaComplex[] zDot)
            throws MaxCountExceededException, DimensionMismatchException {
        TPSOperator rho = new TPSOperator(z);
        
        TPSOperator rhodot = new TPSOperator(zDot);
        rhodot.set(0);
        
        // Individual Hamiltonians
        for(int a = 0; a < n; ++a) {
            t1.set(0,0.5*d[a]);
            t2.set(0,-0.5*d[a]);
            rhodot.pauliLeft(PauliOp.Z, a, t1, rho);
            rhodot.pauliRight(PauliOp.Z, a, t2, rho);
        }
        
        // Interaction potentials
        for(int a = 0; a < n; ++a) {
            sum1.set(0);
            for(int b = 0; b < n; ++b) {
                if(a == b) continue;
                t1.set(0,-0.5*gamma*alpha.getImaginary());
                sum1.pauliLeft(PauliOp.MINUS, b, t1, rho);
            }
            t1.set(1,0);
            rhodot.pauliLeft(PauliOp.PLUS, a, t1, sum1);
        }
        
        for(int b = 0; b < n; ++b) {
            sum1.set(0);
            for(int a = 0; a < n; ++a) {
                if(a == b) continue;
                t1.set(0,0.5*gamma*alpha.getImaginary());
                sum1.pauliRight(PauliOp.PLUS, a, t1, rho);
            }
            t1.set(1,0);
            rhodot.pauliRight(PauliOp.MINUS, b, t1, sum1);
        }
        
        // f-terms
        for(int a = 0; a < n; ++a) {
            sum1.set(0);
            for(int b = 0; b < n; ++b) {
                double f = alpha.getReal();
                if(a == b && f != 0) {
                    f = 1.0;
                }
                t1.set(-0.5*gamma*f,0);
                sum1.pauliLeft(PauliOp.MINUS, b, t1, rho);
            }
            t1.set(1,0);
            rhodot.pauliLeft(PauliOp.PLUS, a, t1, sum1);
        }
        
        for(int b = 0; b < n; ++b) {
            sum1.set(0);
            for(int a = 0; a < n; ++a) {
                double f = alpha.getReal();
                if(a == b && f != 0) {
                    f = 1.0;
                }
                t1.set(-0.5*gamma*f,0);
                sum1.pauliRight(PauliOp.PLUS, a, t1, rho);
            }
            t1.set(1,0);
            rhodot.pauliRight(PauliOp.MINUS, b, t1, sum1);
        }
        
        for(int b = 0; b < n; ++b) {
            sum1.set(0);   
            for(int a = 0; a < n; ++a) {
                double f = alpha.getReal();
                if(a == b && f != 0) {
                    f = 1.0;
                }
                t1.set(gamma*f,0);
                sum1.pauliRight(PauliOp.PLUS, a, t1, rho);
            }
            t1.set(1,0);
            rhodot.pauliLeft(PauliOp.MINUS, b, t1, sum1);
        }
        
        // Pumping terms
        for(int a = 0; a < n; ++a) {
            t1.set(-0.5*w,0);
            
            sum1.set(rho);
            sum1.pauliLeft(PauliOp.PLUS, a);
            rhodot.pauliLeft(PauliOp.MINUS, a, t1, sum1);
            
            sum1.set(rho);
            sum1.pauliRight(PauliOp.MINUS, a);
            rhodot.pauliRight(PauliOp.PLUS, a, t1, sum1);
            
            t1.set(w,0);
            rhodot.pauliLeft(PauliOp.PLUS, a, t1, sum1);
        }
    }

    @Override
    public int getDimension() {
        return SynchUtils.pow(2, 2*n);
    }

}
