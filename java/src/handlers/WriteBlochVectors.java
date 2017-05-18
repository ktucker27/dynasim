package handlers;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;

import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.sampling.StepHandler;
import org.apache.commons.math3.ode.sampling.StepInterpolator;

import utils.SynchUtils;

public class WriteBlochVectors implements StepHandler {
    PrintWriter myWriter;
    int[] startIdxArray;
    int n;
    
    public WriteBlochVectors(String filename, int n) throws FileNotFoundException, UnsupportedEncodingException {
        this.n = n;
        startIdxArray = new int[6];
        SynchUtils.getStartIdx(startIdxArray, n);
        myWriter = new PrintWriter(filename, "UTF-8");
    }
    
    @Override
    public void init(double t0, double[] y0, double t) {
        myWriter.print("0, ");
        for(int i = 0; i < n; ++i) {
            myWriter.print(2*y0[2*i] + ", " + 2*y0[2*i+1] + ", " + y0[2*(i+startIdxArray[2])]);
            if(i < n - 1) {
                myWriter.print(", ");
            }
        }
        myWriter.print("\n");
    }
    
    @Override
    public void handleStep(StepInterpolator interpolator, boolean isLast)
            throws MaxCountExceededException {
        double[] y = interpolator.getInterpolatedState();
        
        myWriter.print(interpolator.getInterpolatedTime() + ", ");
        for(int i = 0; i < n; ++i) {
            myWriter.print(2*y[2*i] + ", " + 2*y[2*i+1] + ", " + y[2*(i+startIdxArray[2])]);
            if(i < n - 1) {
                myWriter.print(", ");
            }
        }
        myWriter.print("\n");
        
        if(!isLast) {
            myWriter.flush();
        } else {
            myWriter.close();
        }
    }
}
