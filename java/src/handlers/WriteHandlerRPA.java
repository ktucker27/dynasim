package handlers;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;

import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.sampling.StepHandler;
import org.apache.commons.math3.ode.sampling.StepInterpolator;

import eval.RPAEval;

public class WriteHandlerRPA implements StepHandler {
    PrintWriter myWriter;
    double myTmin;
    int n;
    
    public WriteHandlerRPA(String filename, int n) throws FileNotFoundException, UnsupportedEncodingException {
        myWriter = new PrintWriter(filename, "UTF-8");
        myTmin = 0.0;
        this.n = n;
    }
    
    @Override
    public void init(double t0, double[] y0, double t) {
        if(myTmin == 0.0) {
            RPAEval eval = new RPAEval(n);
            eval.setVals(y0);
            
            myWriter.print("0, ");
            myWriter.print(eval.compCorr());
            myWriter.print("\n");
        }
    }

    @Override
    public void handleStep(StepInterpolator interpolator, boolean isLast) throws MaxCountExceededException {
        double[] y = interpolator.getInterpolatedState();
        double t = interpolator.getInterpolatedTime();
        
        if(t >= myTmin) {
            RPAEval eval = new RPAEval(n);
            eval.setVals(y);
            
            myWriter.print(t + ", ");
            myWriter.print(eval.compCorr());
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
