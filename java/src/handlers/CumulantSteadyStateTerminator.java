package handlers;

import java.util.ArrayList;

import org.apache.commons.math3.ode.events.EventHandler;
import org.apache.commons.math3.ode.sampling.StepHandler;
import org.apache.commons.math3.ode.sampling.StepInterpolator;

import utils.SynchUtils;

public class CumulantSteadyStateTerminator implements EventHandler {

    int n;
    double myStopTime;
    SteadyStateDetector myDetector;
    int myMaxIterations;
    boolean mySteadyStateReached;
    boolean myQuietMode;
    
    private class SteadyStateDetector implements StepHandler {
        double myMinTime;
        double myTimeDelta;
        double myElapsedTime;
        double myTol;
        int myNumSteps;
        int myTotalSteps;
        ArrayList<Double> mySpinCorr;
        
        public SteadyStateDetector(double minTime, double timeDelta, double tol, int numSteps) {
            myMinTime = minTime;
            myTimeDelta = timeDelta;
            myElapsedTime = 0.0;
            myTol = tol;
            myNumSteps = numSteps;
            myTotalSteps = 0;
            mySpinCorr = new ArrayList<Double>(numSteps);
        }
        
        @Override
        public void init(double t0, double[] y0, double t) {
            mySpinCorr.add(SynchUtils.compCorr(y0, n).getReal());
        }

        @Override
        public void handleStep(StepInterpolator interpolator, boolean isLast) {
            ++myTotalSteps;
            
            if(myStopTime >= 0) return;

            if(myTotalSteps >= myMaxIterations) {
                myStopTime = interpolator.getInterpolatedTime();
                System.err.println("Maximum number of iterations (" + myMaxIterations + ") exceeded at time " + myStopTime);
                return;
            }

            if(interpolator.getInterpolatedTime() - myElapsedTime < myTimeDelta - 1.0e-5) {
                return;
            }
            
            myElapsedTime = interpolator.getInterpolatedTime();
            
            // Copy the new state into the solution stack
            if(mySpinCorr.size() >= myNumSteps) {
                mySpinCorr.remove(0);
            }
            
            double[] y = interpolator.getInterpolatedState(); 
            mySpinCorr.add(SynchUtils.compCorr(y, n).getReal());
            
            // Determine if we have achieved steady state
            if(mySpinCorr.size() >= myNumSteps && 
               interpolator.getInterpolatedTime() >= myMinTime &&
               SynchUtils.isSteadyState(mySpinCorr, myNumSteps, myTol)) {
                myStopTime = interpolator.getInterpolatedTime();
                mySteadyStateReached = true;
                if(!myQuietMode) System.out.println("Steady state detected at time " + myStopTime);
            }
        }
    }
    
    public CumulantSteadyStateTerminator(double minTime, double timeDelta, int numSteps, int maxIterations, double tol, int n) {
        this.n = n;
        myStopTime = -1.0;
        myMaxIterations = maxIterations;
        mySteadyStateReached = false;
        myDetector = new SteadyStateDetector(minTime, timeDelta, tol, numSteps);
        myQuietMode = false;
    }
    
    public void setQuietMode(boolean quiet) {
        myQuietMode = quiet;
    }
    
    public StepHandler getDetector() {
        return myDetector;
    }
    
    @Override
    public void init(double t0, double[] y0, double t) {
    }

    @Override
    public double g(double t, double[] y) {
        if(myStopTime < 0) return 1.0;
        
        double g = 1.0;
        if(t > myStopTime) {
            g = -1.0;
        } else if(t == myStopTime) {
            g = 0.0;
        }
        
        return g;
    }

    @Override
    public Action eventOccurred(double t, double[] y, boolean increasing) {
        if(!myQuietMode) System.out.println("STOP EVENT");
        return EventHandler.Action.STOP;
    }

    @Override
    public void resetState(double t, double[] y) {
    }
    
    public boolean getSteadyStateReached() {
        return mySteadyStateReached;
    }
    
    public double getStopTime() {
        return myStopTime;
    }

}
