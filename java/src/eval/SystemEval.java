package eval;


/**
 * 
 * @author kristophertucker
 *
 * Contains the logic for transforming an array of doubles into various values of physical interest.
 * Evaluators can also provide values needed for specific dynamical systems, as well as initial
 * conditions
 * 
 */
public interface SystemEval {
    
    public enum InitAngleType {
        EQUAL_SPACING,
        RANDOM,
        CONST
    }

    public int getN();
    
    public int getRealDimension();

    public int getDimension();
    
    public double getOrderParam(double[] y);
    
    public double getAvgSigmaz(double[] y);

    public double getAvgSigmazz(double[] y);
    
    public void getBlochVectors(double[] y, double[] xs, double[] ys, double[] zs);
    
    public void initSpinUpX(double[] y0);
    
    public void initialize(double[] y0, double zenith, double phase, InitAngleType zenithType, InitAngleType phaseType);
}
