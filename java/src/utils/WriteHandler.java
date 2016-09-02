package utils;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;

import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.sampling.StepHandler;
import org.apache.commons.math3.ode.sampling.StepInterpolator;

public class WriteHandler implements StepHandler {
    PrintWriter myWriter;
    int myWriteIdx;
    
    public WriteHandler(String filename, int writeIdx) throws FileNotFoundException, UnsupportedEncodingException {
        myWriteIdx = writeIdx;
        myWriter = new PrintWriter(filename, "UTF-8");
    }
    
    @Override
    public void init(double t0, double[] y0, double t) {
        myWriter.print("0, " + y0[myWriteIdx] + "\n");
    }
    
    @Override
    public void handleStep(StepInterpolator interpolator, boolean isLast)
            throws MaxCountExceededException {
        double[] y = interpolator.getInterpolatedState();
        
        myWriter.print(interpolator.getInterpolatedTime() + ", ");
        myWriter.print(y[myWriteIdx]);
        myWriter.print("\n");
        
        if(!isLast) {
            myWriter.flush();
        } else {
            myWriter.close();
        }
    }
}
