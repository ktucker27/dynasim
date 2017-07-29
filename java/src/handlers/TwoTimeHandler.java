package handlers;

import java.io.ByteArrayOutputStream;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;

import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.sampling.StepHandler;
import org.apache.commons.math3.ode.sampling.StepInterpolator;

import utils.DynaComplex;

public class TwoTimeHandler implements StepHandler {
    PrintWriter myWriter;
    ByteArrayOutputStream myBuffer;
    String myFilename;
    DynaComplex mySum;
    DynaComplex t1;
    
    public TwoTimeHandler(String filename, boolean liveStream) throws FileNotFoundException, UnsupportedEncodingException {
        mySum = new DynaComplex();
        t1 = new DynaComplex();
        myFilename = filename;
        
        if(liveStream) {
            myBuffer = null;
            myWriter = new PrintWriter(filename, "UTF-8");
        } else {
            myBuffer = new ByteArrayOutputStream();
            myWriter = new PrintWriter(myBuffer);
        }
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
        
        if(!isLast) {
            myWriter.flush();
        } else {
            myWriter.close();
            if(myBuffer != null) {
                try {
                    PrintWriter fileWriter = new PrintWriter(myFilename, "UTF-8");
                    fileWriter.print(myBuffer.toString());
                    fileWriter.close();
                } catch(FileNotFoundException ex) {
                    System.err.println("Could not open output file " + myFilename);
                    System.err.println(ex.getMessage());
                } catch(UnsupportedEncodingException ex) {
                    System.err.println(ex.getMessage());
                }
            }
        }
    }

    private void compBigOlSum(double[] y) {
        int n = y.length;
        mySum.set(0, 0);
        for(int i = 0; i < n; i += 2) {
            mySum.add(t1.set(y[i], y[i+1]));
        }
    }
}
