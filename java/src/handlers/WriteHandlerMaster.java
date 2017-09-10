package handlers;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;

import ode.DynaComplexODEAdapter;

import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.sampling.StepHandler;
import org.apache.commons.math3.ode.sampling.StepInterpolator;

import utils.DynaComplex;
import utils.SynchUtils;
import utils.TPSOperator;
import utils.TPSOperator.PauliOp;

public class WriteHandlerMaster implements StepHandler {
    PrintWriter myWriter;
    TPSOperator tmp;
    DynaComplex t1;
    DynaComplex[] z;
    double myTmin;
    int n;
    
    public WriteHandlerMaster(String filename, int n) throws FileNotFoundException, UnsupportedEncodingException {
        myWriter = new PrintWriter(filename, "UTF-8");
        myTmin = 0.0;
        this.n = n;
        tmp = new TPSOperator(n);
        t1 = new DynaComplex(1,0);
        int dim = SynchUtils.pow(2, 2*n);
        z = new DynaComplex[dim];
        for(int i = 0; i < dim; ++i) {
            z[i] = new DynaComplex(0,0);
        }
    }
    
    @Override
    public void init(double t0, double[] y0, double t) {
        if(myTmin == 0.0) {
            myWriter.print("0, ");
            
            DynaComplexODEAdapter.toComplex(y0, z);
            TPSOperator rho = new TPSOperator(z);
            
            tmp.set(rho);
            tmp.pauliLeft(PauliOp.MINUS, 1);
            tmp.pauliLeft(PauliOp.PLUS, 0);
            tmp.trace(t1);
            myWriter.print(t1.getReal() + ", ");
            myWriter.print(t1.getImaginary() + ", ");
            
            t1.set(1,0);
            tmp.set(0);
            tmp.pauliRight(PauliOp.PLUS, 0, t1, rho);
            tmp.trace(t1);
            myWriter.print(t1.getReal() + ", ");
            myWriter.print(t1.getImaginary() + ", ");
            
            t1.set(1,0);
            tmp.set(0);
            tmp.pauliRight(PauliOp.Z, 0, t1, rho);
            tmp.trace(t1);
            myWriter.print(t1.getReal() + ", ");
            myWriter.print(t1.getImaginary());
            
            myWriter.print("\n");
        }
    }

    @Override
    public void handleStep(StepInterpolator interpolator, boolean isLast)
            throws MaxCountExceededException {
        double[] y = interpolator.getInterpolatedState();
        double t = interpolator.getInterpolatedTime();
        
        if(t >= myTmin) {
            myWriter.print(t + ", ");
            
            DynaComplexODEAdapter.toComplex(y, z);
            TPSOperator rho = new TPSOperator(z);
            
            tmp.set(rho);
            tmp.pauliLeft(PauliOp.MINUS, 1);
            tmp.pauliLeft(PauliOp.PLUS, 0);
            tmp.trace(t1);
            myWriter.print(t1.getReal() + ", ");
            myWriter.print(t1.getImaginary() + ", ");
            
            t1.set(1,0);
            tmp.set(0);
            tmp.pauliRight(PauliOp.PLUS, 0, t1, rho);
            tmp.trace(t1);
            myWriter.print(t1.getReal() + ", ");
            myWriter.print(t1.getImaginary() + ", ");
            
            t1.set(1,0);
            tmp.set(0);
            tmp.pauliRight(PauliOp.Z, 0, t1, rho);
            tmp.trace(t1);
            myWriter.print(t1.getReal() + ", ");
            myWriter.print(t1.getImaginary());
            
            myWriter.print("\n");

            if(!isLast) {
                myWriter.flush();
            }
        }
        
        if(isLast) {
            myWriter.close();
        }
    }
    
    public void setMinTime(double tmin) {
        myTmin = tmin;
    }
}
