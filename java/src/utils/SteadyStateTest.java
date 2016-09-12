package utils;

import java.util.ArrayDeque;

public interface SteadyStateTest {
    public void init(double t0, double[] y0, double t);
    
    public boolean isSteadyState(ArrayDeque<double[]> soln, int totalSteps);
}
