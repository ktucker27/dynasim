package eval;

import ode.DynaComplexODEAdapter;
import utils.DynaComplex;
import utils.SynchUtils;
import utils.TPSOperator;
import utils.TPSOperator.PauliOp;

public class MasterEval implements SystemEval {

    int n;
    TPSOperator myRho;
    DynaComplex[] z;
    DynaComplex t1;
    
    public MasterEval(int n) {
        this.n = n;
        myRho = new TPSOperator(n);
        t1 = new DynaComplex();
        
        z = new DynaComplex[getDimension()];
        for(int i = 0; i < z.length; ++i) {
            z[i] = new DynaComplex();
        }
    }
    
    @Override
    public int getN() {
        return n;
    }

    @Override
    public int getRealDimension() {
        return 2*getDimension();
    }
    
    @Override
    public int getDimension() {
        return SynchUtils.pow(2, 2*n);
    }

    @Override
    public double getOrderParam(double[] y) {
        DynaComplexODEAdapter.toComplex(y, z);
        TPSOperator rho = new TPSOperator(z);
        
        myRho.set(rho);
        myRho.pauliLeft(PauliOp.MINUS, 1);
        myRho.pauliLeft(PauliOp.PLUS, 0);
        myRho.trace(t1);
        
        return t1.getReal();
    }

    @Override
    public void initSpinUpX(double[] y0) {
        myRho.set(1.0/Math.pow(2,n));
        DynaComplexODEAdapter.toReal(myRho.getVals(), y0);
    }
}
