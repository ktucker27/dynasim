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
        DynaComplex sigmap = SynchUtils.compSigmapAvg(y0, n);
        myWriter.print(SynchUtils.compCorr(y0, n).getReal() + ", ");
        myWriter.print(sigmap.getReal() + ", ");
        myWriter.print(sigmap.getImaginary() + ", ");
        myWriter.print(SynchUtils.compSigmazAvg(y0, n));
        myWriter.print("\n");
    }

    @Override
    public void handleStep(StepInterpolator interpolator, boolean isLast)
            throws MaxCountExceededException {
        double[] y = interpolator.getInterpolatedState();
        
        DynaComplex sigmap = SynchUtils.compSigmapAvg(y, n);
        myWriter.print(interpolator.getInterpolatedTime() + ", ");
        myWriter.print(SynchUtils.compCorr(y, n).getReal() + ", ");
        myWriter.print(sigmap.getReal() + ", ");
        myWriter.print(sigmap.getImaginary() + ", ");
        myWriter.print(SynchUtils.compSigmazAvg(y, n));
        myWriter.print("\n");
        
        if(!isLast) {
            myWriter.flush();
        } else {
            myWriter.close();
        }
    }
}
