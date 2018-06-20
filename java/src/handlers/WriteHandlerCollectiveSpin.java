package handlers;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;

import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.sampling.StepHandler;
import org.apache.commons.math3.ode.sampling.StepInterpolator;

import eval.SystemEval;
import utils.DynaComplex;

public class WriteHandlerCollectiveSpin implements StepHandler {
    PrintWriter myWriter;
    SystemEval myEval;
    double[] myEs;
    DynaComplex[][] myEss;
    
    public WriteHandlerCollectiveSpin(String filename, SystemEval eval) throws FileNotFoundException, UnsupportedEncodingException {
        myWriter = new PrintWriter(filename, "UTF-8");
        myEval = eval;
        
        myEs = new double[3];
        
        myEss = new DynaComplex[3][3];
        for(int i = 0; i < 3; ++i) {
            for(int j = 0; j < 3; ++j) {
                myEss[i][j] = new DynaComplex(0,0);
            }
        }
    }
    
    @Override
    public void init(double t0, double[] y0, double t) {        
        printVals(t0, y0);
    }

    @Override
    public void handleStep(StepInterpolator interpolator, boolean isLast)
            throws MaxCountExceededException {
        double[] y = interpolator.getInterpolatedState();
        double t = interpolator.getInterpolatedTime();
        
        printVals(t, y);
        
        if(!isLast) {
            myWriter.flush();
        } else {
            myWriter.close();
        }
    }
    
    private void printVals(double t, double[] y) {
        myEval.getFirstOrderCollectiveEvs(y, myEs);
        myEval.getSecondOrderCollectiveEvs(y, myEss);
        
        myWriter.print(t);
        for(int i = 0; i < 3; ++i) {
            myWriter.print(", " + myEs[i]);
        }
        
        for(int i = 0; i < 3; ++i) {
            for(int j = 0; j < 3; ++j) {
                myWriter.print(", " + myEss[i][j].getReal());
                myWriter.print(", " + myEss[i][j].getImaginary());
            }
        }
        myWriter.print("\n");
    }
}
