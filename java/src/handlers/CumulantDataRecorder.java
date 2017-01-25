package handlers;

import java.util.ArrayList;

import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.sampling.StepHandler;
import org.apache.commons.math3.ode.sampling.StepInterpolator;

import utils.SynchUtils;

public class CumulantDataRecorder implements StepHandler {

    int n;
    double myTimeStep;
    double myElapsedTime;
    double myTol;
    ArrayList<Double> myTimes;
    ArrayList<Double> mySpinCorr;
    
    public CumulantDataRecorder(double timeStep, double totalTime, int n) {
        this.n = n;
        myTimeStep = timeStep;
        myElapsedTime = 0.0;
        myTol = 1.0e-5;
        int estSize = (int) (Math.ceil(totalTime/timeStep) + 1.5);
        myTimes = new ArrayList<Double>(estSize);
        mySpinCorr = new ArrayList<Double>(estSize);
    }
    
    @Override
    public void init(double t0, double[] y0, double t) {
        myTimes.add(0.0);
        mySpinCorr.add(SynchUtils.compCorr(y0, n).getReal());
    }

    @Override
    public void handleStep(StepInterpolator interpolator, boolean isLast)
            throws MaxCountExceededException {
        double time = interpolator.getInterpolatedTime();
        if(time - myElapsedTime >= myTimeStep - myTol) {
            double[] y = interpolator.getInterpolatedState();
            myTimes.add(time);
            mySpinCorr.add(SynchUtils.compCorr(y, n).getReal());
            myElapsedTime = time;
        }
    }
    
    public ArrayList<Double> getTimes() {
        return myTimes;
    }
    
    public ArrayList<Double> getSpinCorr() {
        return mySpinCorr;
    }
}
