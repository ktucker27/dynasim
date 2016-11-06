package coupling;

import org.apache.commons.math3.complex.Complex;

public class ConstCoupling implements ComplexCoupling {
    
    double f;
    double g;
    Complex alpha;
    
    public ConstCoupling(double f, double g) {
        this.f = f;
        this.g = g;
        alpha = new Complex(f, g);
    }
    
    @Override
    public double getF(int row, int col) {
        return f;
    }

    @Override
    public double getG(int row, int col) {
        return g;
    }

    @Override
    public Complex getAlpha(int row, int col) {
        return alpha;
    }

}
