package handlers;

import java.util.ArrayDeque;
import java.util.Iterator;

import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;

import utils.OrderParameterSolution;

public class MeanFieldSteadyStateTest implements SteadyStateTest {
    
    double myTol;
    double mySum;
    
    public MeanFieldSteadyStateTest(double tol) {
        myTol = tol;
        mySum = 0.0;
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
            
            double z = OrderParameterSolution.compOrderParam(currState);
            zvals[idx] = z;
            mySum += z;
            ++idx;
        }
        
        double mean = mySum/(double)totalSteps;
        
        if(stdDev.evaluate(zvals)/mean >= myTol) {
            isSteadyState = false;
        }
        
        return isSteadyState;
    }

}
