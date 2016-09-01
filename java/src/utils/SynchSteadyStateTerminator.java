package utils;

import java.util.ArrayDeque;
import java.util.Iterator;

import org.apache.commons.math3.ode.events.EventHandler;
import org.apache.commons.math3.ode.sampling.StepHandler;
import org.apache.commons.math3.ode.sampling.StepInterpolator;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;

public class SynchSteadyStateTerminator implements EventHandler {

    double myStopTime;
    SteadyStateDetector myDetector;
    int myMaxIterations;
    boolean myQuietMode;
    
    private class SteadyStateDetector implements StepHandler {
        int myNumSteps;
        int myTotalSteps;
        ArrayDeque<double[]> mySoln;
        int n;
        int myMargin;
        
        public SteadyStateDetector(int numSteps) {
            myNumSteps = numSteps;
            myTotalSteps = 0;
        }
        
        @Override
        public void init(double t0, double[] y0, double t) {
            mySoln = new ArrayDeque<double[]>(myNumSteps);
            n = y0.length/3;
            myMargin = (int) (0.1*n) + 1;
        }

        @Override
        public void handleStep(StepInterpolator interpolator, boolean isLast) {
            ++myTotalSteps;
            
            if(myStopTime >= 0) return;

            if(myTotalSteps >= myMaxIterations) {
                myStopTime = interpolator.getInterpolatedTime();
                System.out.println("Maximum number of iterations (" + myMaxIterations + ") exceeded at time " + myStopTime);
                return;
            }

            // Copy the new state into the solution stack
            if(mySoln.size() >= myNumSteps) {
                mySoln.pollLast();
            }
            
            double[] state = new double[3*n];
            for(int i = 0; i < 3*n; ++i) {
                state[i] = interpolator.getInterpolatedState()[i];
            }
            mySoln.push(state);

            if(myTotalSteps >= myNumSteps) {
                // Determine if we have achieved steady state
                boolean isSteadyState = true;
                StandardDeviation stdDev = new StandardDeviation();
                Iterator<double[]> iter = mySoln.descendingIterator();
                double[] rvals = new double[myNumSteps];
                int idx = 0;
                while(iter.hasNext()) {
                    double[] currState = iter.next();
                    
                    double[] rvec = new double[n - 2*myMargin];
                    for(int i = 0; i < rvec.length; ++i) {
                        rvec[i] = Math.abs(currState[n+i+myMargin]);
                    }
                    
                    if(!myQuietMode && myTotalSteps % 1000 == 0 && idx % 100 == 0) {
                        System.out.println("rvec: " + stdDev.evaluate(rvec));
                    }
                    
                    if(stdDev.evaluate(rvec) > 5.0e-4) {
                        isSteadyState = false;
                        break;
                    }
                    
                    rvals[idx] = currState[n+n/2];
                    ++idx;
                }
                
                if(isSteadyState) {
                    if(!myQuietMode && myTotalSteps % 1000 == 0) {
                        System.out.println("rvals: " + stdDev.evaluate(rvals));
                    }
                    
                    if(stdDev.evaluate(rvals) < 1.0e-4) {
                        myStopTime = interpolator.getInterpolatedTime();
                        System.out.println("Steady state detected at time " + myStopTime);
                    }
                }
            }
        }        
    }
    
    public SynchSteadyStateTerminator(int numSteps, int maxIterations) {
       myStopTime = -1.0;
       myMaxIterations = maxIterations;
       myQuietMode = true;
       myDetector = new SteadyStateDetector(numSteps);
    }
    
    public void setQuietMode(boolean quietMode) {
        myQuietMode = quietMode;
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

}
