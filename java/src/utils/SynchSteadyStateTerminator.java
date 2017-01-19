package utils;

import java.util.ArrayDeque;

import org.apache.commons.math3.ode.events.EventHandler;
import org.apache.commons.math3.ode.sampling.StepHandler;
import org.apache.commons.math3.ode.sampling.StepInterpolator;

public class SynchSteadyStateTerminator implements EventHandler {

    double myStopTime;
    SteadyStateDetector myDetector;
    int myMaxIterations;
    boolean myReturnOk;
    
    private class SteadyStateDetector implements StepHandler {
        double myMinTime;
        int myNumSteps;
        int myTotalSteps;
        ArrayDeque<double[]> mySoln;
        SteadyStateTest myTest;
        
        public SteadyStateDetector(SteadyStateTest test, double minTime, int numSteps) {
            myMinTime = minTime;
            myNumSteps = numSteps;
            myTotalSteps = 0;
            myTest = test;
        }
        
        @Override
        public void init(double t0, double[] y0, double t) {
            mySoln = new ArrayDeque<double[]>(myNumSteps);
            myTest.init(t0, y0, t);
        }

        @Override
        public void handleStep(StepInterpolator interpolator, boolean isLast) {
            ++myTotalSteps;
            
            if(myStopTime >= 0) return;

            if(myTotalSteps >= myMaxIterations) {
                myStopTime = interpolator.getInterpolatedTime();
                System.out.println("Maximum number of iterations (" + myMaxIterations + ") exceeded at time " + myStopTime);
                myReturnOk = false;
                return;
            }

            // Copy the new state into the solution stack
            if(mySoln.size() >= myNumSteps) {
                mySoln.pollLast();
            }
            
            int n = interpolator.getInterpolatedState().length;
            double[] state = new double[n];
            for(int i = 0; i < n; ++i) {
                state[i] = interpolator.getInterpolatedState()[i];
            }
            mySoln.push(state);

            // Determine if we have achieved steady state
            if(myTotalSteps >= myNumSteps && 
               myTest.isSteadyState(mySoln, myTotalSteps) &&
               interpolator.getInterpolatedTime() >= myMinTime) {
                myStopTime = interpolator.getInterpolatedTime();
                System.out.println("Steady state detected at time " + myStopTime);
            }
        }        
    }
    
    public SynchSteadyStateTerminator(SteadyStateTest test, double minTime, int numSteps, int maxIterations) {
       myStopTime = -1.0;
       myMaxIterations = maxIterations;
       myReturnOk = true;
       myDetector = new SteadyStateDetector(test, minTime, numSteps);
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
        System.out.println("STOP EVENT");
        return EventHandler.Action.STOP;
    }

    @Override
    public void resetState(double t, double[] y) {
    }
    
    public boolean getReturnOk() {
        return myReturnOk;
    }

}
