package handlers;

import java.util.ArrayDeque;
import java.util.Deque;
import java.util.Iterator;

import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.sampling.StepHandler;
import org.apache.commons.math3.ode.sampling.StepInterpolator;

import eval.SystemEval;

public class DataRecorder implements StepHandler {

    SystemEval myEval;
    double myTimeStep;
    double myElapsedTime;
    double myTol;
    int myMaxSize;
    Deque<Double> myTimes;
    Deque<Double> myOrderParams;
    Deque<Double> myAvgZs;
    
    public DataRecorder(SystemEval eval, double timeStep, double totalTime) {
        myEval = eval;
        myTimeStep = timeStep;
        myElapsedTime = 0.0;
        myTol = 1.0e-5;
        myMaxSize = (int) (Math.ceil(totalTime/timeStep) + 1.5);
        myTimes = new ArrayDeque<Double>(myMaxSize);
        myOrderParams = new ArrayDeque<Double>(myMaxSize);
        myAvgZs = new ArrayDeque<Double>(myMaxSize);
    }
    
    @Override
    public void init(double t0, double[] y0, double t) {
        myTimes.add(0.0);
        myOrderParams.add(myEval.getOrderParam(y0));
        myAvgZs.add(myEval.getAvgSigmaz(y0));
    }

    @Override
    public void handleStep(StepInterpolator interpolator, boolean isLast)
            throws MaxCountExceededException {
        double time = interpolator.getInterpolatedTime();
        if(time - myElapsedTime >= myTimeStep - myTol) {
            double[] y = interpolator.getInterpolatedState();
            if(myTimes.size() == myMaxSize) {
                dumpStep();
            }
            myTimes.add(time);
            myOrderParams.add(myEval.getOrderParam(y));
            myAvgZs.add(myEval.getAvgSigmaz(y));
            myElapsedTime = time;
        }
    }
    
    public Deque<Double> getOrderParams() {
        return myOrderParams;
    }
    
    public double getMeanOrderParam() {
        return getDequeMean(myOrderParams);
    }

    public Deque<Double> getAvgZs() {
        return myAvgZs;
    }
    
    public double getMeanAvgZs() {
        return getDequeMean(myAvgZs);
    }
    
    private double getDequeMean(Deque<Double> deque) {
        double sum = 0.0;
        Iterator<Double> iter = deque.iterator();
        while(iter.hasNext()) {
            sum += iter.next();
        }
        
        return sum/(double)deque.size();
    }

    private void dumpStep() {
        myTimes.removeFirst();
        myOrderParams.removeFirst();
        myAvgZs.removeFirst();
    }
}
