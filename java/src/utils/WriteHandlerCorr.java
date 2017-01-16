package utils;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;

import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.sampling.StepHandler;
import org.apache.commons.math3.ode.sampling.StepInterpolator;

public class WriteHandlerCorr implements StepHandler {
    PrintWriter myWriter;
    int n;
    
    public WriteHandlerCorr(String filename, int n) throws FileNotFoundException, UnsupportedEncodingException {
        myWriter = new PrintWriter(filename, "UTF-8");
        this.n = n;
    }
    
    @Override
    public void init(double t0, double[] y0, double t) {
        myWriter.print("0, ");
        DynaComplex sigmap = compSigmapAvg(y0, n);
        myWriter.print(compCorr(y0, n).getReal() + ", ");
        myWriter.print(sigmap.getReal() + ", ");
        myWriter.print(sigmap.getImaginary() + ", ");
        myWriter.print(compSigmazAvg(y0, n));
        myWriter.print("\n");
    }

    @Override
    public void handleStep(StepInterpolator interpolator, boolean isLast)
            throws MaxCountExceededException {
        double[] y = interpolator.getInterpolatedState();
        
        DynaComplex sigmap = compSigmapAvg(y, n);
        myWriter.print(interpolator.getInterpolatedTime() + ", ");
        myWriter.print(compCorr(y, n).getReal() + ", ");
        myWriter.print(sigmap.getReal() + ", ");
        myWriter.print(sigmap.getImaginary() + ", ");
        myWriter.print(compSigmazAvg(y, n));
        myWriter.print("\n");
        
        if(!isLast) {
            myWriter.flush();
        } else {
            myWriter.close();
        }
    }

    public static DynaComplex compCorr(double[] y0, int n) {
        double sum = 0.0;
        int startIdx = 2*(2*n + n*(n-1));
        
        for(int i = startIdx; i < startIdx + n*(n-1); i += 2) {
            sum += 2.0*y0[i];
        }
        sum *= 1.0/(n*(n-1));
        
        return new DynaComplex(sum);
    }
    
    public static double compSigmazAvg(double[] y0, int n) {
        double sum = 0.0;
        int startIdx = 2*n*n;
        
        for(int i = startIdx; i < startIdx + 2*n; i += 2) {
            sum += y0[i];
        }
        sum *= 1.0/n;

        return sum;
    }
    
    public static DynaComplex compSigmapAvg(double[] y0, int n) {
        DynaComplex sum = new DynaComplex(0.0);
        DynaComplex z = new DynaComplex();
        int startIdx = 0;
        
        for(int i = startIdx; i < startIdx + 2*n; i += 2) {
            z.set(y0[i], y0[i+1]);
            sum.add(z);
        }
        sum.multiply(1.0/n);

        return sum;
    }
}
