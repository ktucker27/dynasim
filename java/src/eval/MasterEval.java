package eval;

import java.io.FileNotFoundException;
import java.io.PrintWriter;

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
        
        double sum = 0.0;
        for(int a = 0; a < n; ++a) {
            for(int b = a + 1; b < n; ++b) {
                myRho.set(rho);
                myRho.pauliLeft(PauliOp.MINUS, b);
                myRho.pauliLeft(PauliOp.PLUS, a);
                myRho.trace(t1);
                sum += t1.getReal();
            }
        }
        
        return sum/(0.5*n*(n-1));
    }
    
    @Override
    public double getAvgSigmaz(double[] y) {
        DynaComplexODEAdapter.toComplex(y, z);
        TPSOperator rho = new TPSOperator(z);
        
        double sum = 0.0;
        for(int a = 0; a < n; ++a) {
            myRho.set(rho);
            myRho.pauliLeft(PauliOp.Z, a);
            myRho.trace(t1);
            sum += t1.getReal();
        }
        
        return sum/(double)n;
    }
    
    @Override
    public double getAvgSigmazz(double[] y) {
        DynaComplexODEAdapter.toComplex(y, z);
        TPSOperator rho = new TPSOperator(z);
        
        double sum = 0.0;
        for(int a = 0; a < n; ++a) {
            for(int b = a+1; b < n; ++b) {
                myRho.set(rho);
                myRho.pauliLeft(PauliOp.Z, b);
                myRho.pauliLeft(PauliOp.Z, a);
                myRho.trace(t1);
                sum += t1.getReal();
            }
        }
        
        return sum/(double)(n*(n-1)/2);
    }
    
    @Override
    public void getBlochVectors(double[] y, double[] xs, double[] ys, double[] zs) {
        // TODO
        throw new UnsupportedOperationException("Bloch vector retrieval not yet implemented for MasterEval");
    }
    
    @Override
    public void getFirstOrderCollectiveEvs(double[] y, double[] es) {
        throw new UnsupportedOperationException("MasterEval has not implemented getFirstOrderEvs");
    }
    
    @Override
    public void getSecondOrderCollectiveEvs(double[] y, DynaComplex[][] es) {
        throw new UnsupportedOperationException("MasterEval has not implemented getSecondOrderEvs");
    }

    @Override
    public void initSpinUpX(double[] y0) {
        myRho.set(1.0/Math.pow(2,n));
        DynaComplexODEAdapter.toReal(myRho.getVals(), y0);
    }
    
    @Override
    public void initialize(double[] y0, double zenith, double phase, InitAngleType zenithType, InitAngleType phaseType) {
        // TODO
        throw new UnsupportedOperationException("Complete initialization functionality not yet implemented for the master solver");
    }
    
    public void writeZDist(double[] y, String filepath) throws FileNotFoundException {
        DynaComplexODEAdapter.toComplex(y, z);
        TPSOperator rho = new TPSOperator(z);
        
        double[] dist = new double[2*n + 1];
        for(int i = 0; i < dist.length; ++i) {
            dist[i] = 0.0;
        }
        
        for(int i = 0; i < SynchUtils.pow(2, n); ++i) {
            int k = 0;
            int ii = i;
            for(int j = 0; j < n; ++j) {
                k += (ii % 2);
                ii = ii >> 1;
            }
            int eval = n - 2*k;
            dist[eval + n] += rho.getVal(i, i).getReal();
        }
        
        PrintWriter writer = new PrintWriter(filepath);
        for(int i = 0; i < dist.length; ++i) {
            writer.write(dist[i] + "\n");
        }
        writer.close();
    }
    
    public void writeTriples(double[] y) {
        DynaComplexODEAdapter.toComplex(y, z);
        TPSOperator rho = new TPSOperator(z);
        
        // Woefully inefficient for a number of reasons, but will work for
        // the intended experiment
        double[] vals = new double[27];
        int idx = 0;
        for(int al = 0; al < 3; ++al) {
            for(int bt = 0; bt < 3; ++bt) {
                for(int gm = 0; gm < 3; ++gm) {
                    myRho.set(rho);
                    myRho.pauliLeft(PauliOp.values()[gm+3], 0);
                    myRho.pauliLeft(PauliOp.values()[bt+3], 1);
                    myRho.pauliLeft(PauliOp.values()[al+3], 2);
                    myRho.trace(t1);
                    vals[idx] = t1.getReal();
                    if(vals[idx] > 1.0e-10) {
                        System.out.println(al + ", " + bt + ", " + gm + ": " + vals[idx]);
                    }
                }
            }
        }
    }
}
