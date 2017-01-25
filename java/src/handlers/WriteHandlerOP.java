package handlers;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;

import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.sampling.StepHandler;
import org.apache.commons.math3.ode.sampling.StepInterpolator;

import utils.OrderParameterSolution;

public class WriteHandlerOP implements StepHandler {
    PrintWriter myWriter;
    
    public WriteHandlerOP(String filename) throws FileNotFoundException, UnsupportedEncodingException {
        myWriter = new PrintWriter(filename, "UTF-8");
    }
    
    @Override
    public void init(double t0, double[] y0, double t) {
        myWriter.print("0, " + OrderParameterSolution.compOrderParam(y0) + "\n");
    }
    
    @Override
    public void handleStep(StepInterpolator interpolator, boolean isLast)
            throws MaxCountExceededException {
        double[] y = interpolator.getInterpolatedState();
        
        myWriter.print(interpolator.getInterpolatedTime() + ", ");
        myWriter.print(OrderParameterSolution.compOrderParam(y));
        myWriter.print("\n");
        
        if(!isLast) {
            myWriter.flush();
        } else {
            myWriter.close();
        }
    }
}
