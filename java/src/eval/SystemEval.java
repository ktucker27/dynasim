package eval;

/**
 * 
 * @author kristophertucker
 *
 * Wraps an array of doubles and provides various values of common physical interest.
 * Evaluators can also provide values needed for specific dynamical systems, as well
 * as initial conditions
 * 
 */
public interface SystemEval {

    public int getN();
    
    public int getRealDimension();

    public int getDimension();
    
    public double getOrderParam(double[] y);
    
    public double getAvgSigmaz(double[] y);
    
    public void initSpinUpX(double[] y0);
}
