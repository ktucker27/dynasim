package eval;

import java.util.Random;

import org.apache.commons.math3.complex.Complex;

import utils.DynaComplex;

public class MeanFieldEval implements SystemEval {
    
    int n;
    
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
    public double getAvgSigmaz(double[] y) {
        double sum = 0.0;
        
        for(int i = 0; i < n; ++i) {
            sum += y[i];
        }
        sum *= 1.0/n;

        return sum;
    }
    
    @Override
    public double getAvgSigmazz(double[] y) {
        double sum = 0.0;
        
        for(int i = 0; i < n; ++i) {
            for(int j = i+1; j < n; ++j) {
                sum += y[i]*y[j];
            }
        }
        sum *= 1.0/(n*(n-1)/2);

        return sum;
    }
    
    @Override
    public void getBlochVectors(double[] y, double[] xs, double[] ys, double[] zs) {
        Complex z = null;
        for(int i = 0; i < n; ++i) {
            z = Complex.I.multiply(y[2*n+i]).exp().multiply(y[n+i]);
            xs[i] = z.getReal();
            ys[i] = z.getImaginary();
            zs[i] = y[i];
        }
    }
    
    @Override
    public void getFirstOrderCollectiveEvs(double[] y, double[] es) {
        throw new UnsupportedOperationException("MeanFieldEval has not implemented getFirstOrderEvs");
    }
    
    @Override
    public void getSecondOrderCollectiveEvs(double[] y, DynaComplex[][] es) {
        throw new UnsupportedOperationException("MeanFieldEval has not implemented getSecondOrderEvs");
    }

    @Override
    public void initSpinUpX(double[] y0) {
        initialize(y0, 0.5*Math.PI, 0.0, InitAngleType.CONST, InitAngleType.CONST);
    }
    
    @Override
    public void initialize(double[] y0, double zenith, double phase, InitAngleType zenithType, InitAngleType phaseType) {
        Random rand = new Random(5);
        double tip = 0.0;
        double stip = 0.0;
        double ctip = 0.0;
        
        double eps = 1.0e-3;
        
        for(int i = 0; i < n; ++i) {
            
            switch(zenithType) {
            case EQUAL_SPACING:
                tip = eps + i*(Math.PI - 2*eps)/(double)(n-1);
                break;
            case RANDOM:
                tip = rand.nextDouble()*Math.PI;
                break;
            case CONST:
                tip = zenith;
                break;
            }
            
            stip = Math.sin(tip);
            ctip = Math.cos(tip);
            
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
