package handlers;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;

import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.sampling.StepHandler;
import org.apache.commons.math3.ode.sampling.StepInterpolator;

import utils.DynaComplex;

public class TwoTimeHandler implements StepHandler {
    PrintWriter myWriter;
    DynaComplex mySum;
    DynaComplex t1;
    
    public TwoTimeHandler(String filename) throws FileNotFoundException, UnsupportedEncodingException {
        myWriter = new PrintWriter(filename, "UTF-8");
        mySum = new DynaComplex();
        t1 = new DynaComplex();
    }
    
    @Override
    public void init(double t0, double[] y0, double t) {
        myWriter.print("0, ");
        compBigOlSum(y0);
        myWriter.print(mySum.getReal() + ", " + mySum.getImaginary() + "\n");
    }

    @Override
    public void handleStep(StepInterpolator interpolator, boolean isLast)
            throws MaxCountExceededException {
        double[] y = interpolator.getInterpolatedState();
        
        myWriter.print(interpolator.getInterpolatedTime() + ", ");
        compBigOlSum(y);
        myWriter.print(mySum.getReal() + ", " + mySum.getImaginary() + "\n");
    }

    private void compBigOlSum(double[] y) {
        int n = y.length;
        mySum.set(0, 0);
        for(int i = 0; i < n; i += 2) {
            mySum.add(t1.set(y[i], y[i+1]));
        }
    }
}
