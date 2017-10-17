package handlers;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.sampling.StepHandler;
import org.apache.commons.math3.ode.sampling.StepInterpolator;

public class WriteHandlerMeanField implements StepHandler {
    PrintWriter myWriter;
    double myTmin;
    int n;
    
    public WriteHandlerMeanField(String filename, int n) throws FileNotFoundException, UnsupportedEncodingException {
        myWriter = new PrintWriter(filename, "UTF-8");
        myTmin = 0.0;
        this.n = n;
    }
    

    @Override
    public void init(double t0, double[] y0, double t) {
        if(myTmin == 0.0) {
            myWriter.print("0, ");
            myWriter.print(compCorr(y0, n) + ", ");
            myWriter.print(avgSigmax(y0, n) + ", ");
            myWriter.print(avgSigmay(y0, n) + ", ");
            myWriter.print(avgSigmaz(y0, n));
            myWriter.print("\n");
        }
    }

    @Override
    public void handleStep(StepInterpolator interpolator, boolean isLast) throws MaxCountExceededException {
        double[] y = interpolator.getInterpolatedState();
        double t = interpolator.getInterpolatedTime();
        
        if(t >= myTmin) {
            myWriter.print(t + ", ");
            myWriter.print(compCorr(y, n) + ", ");
            myWriter.print(avgSigmax(y, n) + ", ");
            myWriter.print(avgSigmay(y, n) + ", ");
            myWriter.print(avgSigmaz(y, n));
            myWriter.print("\n");

            if(!isLast) {
                myWriter.flush();
            }
        }
        
        if(isLast) {
            myWriter.close();
        }
    }
    
    public static double compCorr(double[] y, int n) {
        // Compute the order parameter (|avg_a(<sigma_a^+>)|)
        Complex z = new Complex(0.0);
        for(int i = 0; i < n; ++i) {
            z = z.add(Complex.I.multiply(y[2*n+i]).exp().multiply(0.5*y[n+i]));
        }
        z = z.multiply(1/(double)n);
        
        return z.abs();
    }
    
    // TODO - Port these to MeanFieldEval and use an instance of that instead
    public static double avgSigmax(double[] y, int n) {
        double sum = 0.0;
        for(int i = 0; i < n; ++i) {
            sum += y[n+i]*Math.cos(y[2*n+i]);
        }
        
        return sum/(double)n;
    }
    
    public static double avgSigmay(double[] y, int n) {
        double sum = 0.0;
        for(int i = 0; i < n; ++i) {
            sum += y[n+i]*Math.sin(y[2*n+i]);
        }
        
        return sum/(double)n;
    }
    
    public static double avgSigmaz(double[] y, int n) {
        double sum = 0.0;
        for(int i = 0; i < n; ++i) {
            sum += y[i];
        }
        
        return sum/(double)n;
    }
}
