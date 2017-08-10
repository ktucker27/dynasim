package handlers;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;

import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.sampling.StepHandler;
import org.apache.commons.math3.ode.sampling.StepInterpolator;

import utils.DynaComplex;
import utils.SynchUtils;

public class WriteHandlerCorr implements StepHandler {
    PrintWriter myWriter;
    double myTmin;
    int n;
    
    public WriteHandlerCorr(String filename, int n) throws FileNotFoundException, UnsupportedEncodingException {
        myWriter = new PrintWriter(filename, "UTF-8");
        myTmin = 0.0;
        this.n = n;
    }
    
    @Override
    public void init(double t0, double[] y0, double t) {
        if(myTmin == 0.0) {
            myWriter.print("0, ");
            DynaComplex sigmap = SynchUtils.compSigmapAvg(y0, n);
            DynaComplex sigmazp = SynchUtils.compSigmazpAvg(y0, n);
            DynaComplex sigmapp = SynchUtils.compSigmappAvg(y0, n);
            myWriter.print(SynchUtils.compCorr(y0, n).getReal() + ", ");
            myWriter.print(sigmap.getReal() + ", ");
            myWriter.print(sigmap.getImaginary() + ", ");
            myWriter.print(SynchUtils.compSigmazAvg(y0, n) + ", ");
            myWriter.print(sigmazp.getReal() + ", ");
            myWriter.print(sigmazp.getImaginary() + ", ");
            myWriter.print(SynchUtils.compSigmazzAvg(y0, n) + ", ");
            myWriter.print(sigmapp.getReal() + ", ");
            myWriter.print(sigmapp.getImaginary() + ", ");
            myWriter.print("\n");
        }
    }

    @Override
    public void handleStep(StepInterpolator interpolator, boolean isLast)
            throws MaxCountExceededException {
        double[] y = interpolator.getInterpolatedState();
        double t = interpolator.getInterpolatedTime();
        
        if(t >= myTmin) {
            DynaComplex sigmap = SynchUtils.compSigmapAvg(y, n);
            DynaComplex sigmazp = SynchUtils.compSigmazpAvg(y, n);
            DynaComplex sigmapp = SynchUtils.compSigmappAvg(y, n);
            myWriter.print(t + ", ");
            myWriter.print(SynchUtils.compCorr(y, n).getReal() + ", ");
            myWriter.print(sigmap.getReal() + ", ");
            myWriter.print(sigmap.getImaginary() + ", ");
            myWriter.print(SynchUtils.compSigmazAvg(y, n) + ", ");
            myWriter.print(sigmazp.getReal() + ", ");
            myWriter.print(sigmazp.getImaginary() + ", ");
            myWriter.print(SynchUtils.compSigmazzAvg(y, n) + ", ");
            myWriter.print(sigmapp.getReal() + ", ");
            myWriter.print(sigmapp.getImaginary() + ", ");
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
