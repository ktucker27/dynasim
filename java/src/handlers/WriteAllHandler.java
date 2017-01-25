package handlers;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;

import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.sampling.StepHandler;
import org.apache.commons.math3.ode.sampling.StepInterpolator;

public class WriteAllHandler implements StepHandler {
    PrintWriter myWriter;
    int myDim;
    
    public WriteAllHandler(String filename, int dim) throws FileNotFoundException, UnsupportedEncodingException {
        myDim = dim;
        myWriter = new PrintWriter(filename, "UTF-8");
    }
    
    @Override
    public void init(double t0, double[] y0, double t) {
        myWriter.print("0, ");
        for(int i = 0; i < myDim; ++i) {
            myWriter.print(y0[i]);
            if(i < myDim - 1) {
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
        for(int i = 0; i < myDim; ++i) {
            myWriter.print(y[i]);
            if(i < myDim - 1) {
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
