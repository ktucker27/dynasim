package handlers;

import java.util.ArrayDeque;
import java.util.Iterator;

import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;

import utils.SynchUtils;

public class SpinCorrSteadyStateTest implements SteadyStateTest {
    
    double myTol;
    double mySum;
    int n;
    
    public SpinCorrSteadyStateTest(double tol, int n) {
        myTol = tol;
        mySum = 0.0;
        this.n = n;
    }
    
    @Override
    public void init(double t0, double[] y0, double t) {
    }

    @Override
    public boolean isSteadyState(ArrayDeque<double[]> soln, int totalSteps) {
        int numSteps = soln.size();

        if(totalSteps % numSteps != 0) return false;
        
        boolean isSteadyState = true;
        
        StandardDeviation stdDev = new StandardDeviation();
        Iterator<double[]> iter = soln.descendingIterator();
        double[] zvals = new double[numSteps];
        int idx = 0;
        while(iter.hasNext()) {
            double[] currState = iter.next();
            
            double z = SynchUtils.compCorr(currState, n).getReal();
            zvals[idx] = z;
            mySum += z;
            ++idx;
        }
        
        double mean = mySum/(double)totalSteps;
        
//        System.out.println(stdDev.evaluate(zvals) + " " + mean + " " + stdDev.evaluate(zvals)/mean);
        if(stdDev.evaluate(zvals)/mean >= myTol) {
            isSteadyState = false;
        }
        
        return isSteadyState;
    }

}
