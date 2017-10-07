package eval;

import java.util.Random;

import org.apache.commons.math3.complex.Complex;

public class MeanFieldEval implements SystemEval {
    
    int n;
    
    public enum InitPhaseType {
        EQUAL_SPACING,
        RANDOM,
        CONST
    }
    
    public MeanFieldEval(int n) {
        this.n = n;
    }
    
    @Override
    public int getN() {
        return n;
    }

    @Override
    public int getRealDimension() {
        return getDimension();
    }
    
    @Override
    public int getDimension() {
        return 3*n;
    }

    @Override
    public double getOrderParam(double[] y) {
        // Compute the order parameter (|avg_a(<sigma_a^+>)|)
        Complex z = new Complex(0.0);
        for(int i = 0; i < n; ++i) {
            z = z.add(Complex.I.multiply(y[2*n+i]).exp().multiply(0.5*y[n+i]));
        }
        z = z.multiply(1/(double)n);
        
        return z.abs();
    }

    @Override
    public void initSpinUpX(double[] y0) {
        initialize(y0, 0.5*Math.PI, 0.0, InitPhaseType.CONST);
    }
    
    public void initialize(double[] y0, double tip, double phase, InitPhaseType phaseType) {
        Random rand = new Random(5);
        double stip = Math.sin(tip);
        double ctip = Math.cos(tip);
        
        for(int i = 0; i < n; ++i) {
            
            y0[i] = ctip;
            y0[n+i] = stip;

            switch(phaseType) {
            case EQUAL_SPACING:
                y0[2*n+i] = i*2*Math.PI/(double)n - Math.PI;
                break;
            case RANDOM:
                y0[2*n+i] = rand.nextDouble()*2.0*Math.PI - Math.PI;
                break;
            case CONST:
                y0[2*n+i] = phase;
                break;
            }
        }
    }
}
