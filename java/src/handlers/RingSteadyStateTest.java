package handlers;

import java.util.ArrayDeque;
import java.util.Iterator;

import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;

public class RingSteadyStateTest implements SteadyStateTest {
    int n;
    int myMargin;
    boolean myQuietMode;
    
    public RingSteadyStateTest() {
        n = 0;
        myMargin = 0;
        myQuietMode = true;
    }
    
    public void setQuietMode(boolean quiteMode) {
        myQuietMode = quiteMode;
    }
    
    @Override
    public void init(double t0, double[] y0, double t) {
        n = y0.length/3;
        myMargin = (int) (0.1*n) + 1;
    }
    
    @Override
    public boolean isSteadyState(ArrayDeque<double[]> soln, int totalSteps) {
        boolean isSteadyState = true;
        int numSteps = soln.size();
        StandardDeviation stdDev = new StandardDeviation();
        Iterator<double[]> iter = soln.descendingIterator();
        double[] rvals = new double[numSteps];
        int idx = 0;
        while(iter.hasNext()) {
            double[] currState = iter.next();
            
            double[] rvec = new double[n - 2*myMargin];
            for(int i = 0; i < rvec.length; ++i) {
                rvec[i] = Math.abs(currState[n+i+myMargin]);
            }
            
            if(!myQuietMode && totalSteps % 1000 == 0 && idx % 100 == 0) {
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
            if(!myQuietMode && totalSteps % 1000 == 0) {
                System.out.println("rvals: " + stdDev.evaluate(rvals));
            }
            
            if(stdDev.evaluate(rvals) >= 1.0e-4) {
                isSteadyState = false;
            }
        }
        
        return isSteadyState;
    }
}
