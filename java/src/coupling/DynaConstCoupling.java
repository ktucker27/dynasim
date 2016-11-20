package coupling;

import utils.DynaComplex;

public class DynaConstCoupling implements DynaComplexCoupling {

    double f;
    double g;
    DynaComplex alpha;
    
    public DynaConstCoupling(double f, double g) {
        this.f = f;
        this.g = g;
        alpha = new DynaComplex(f, g);
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
    public DynaComplex getAlpha(int row, int col) {
        return alpha;
    }
}
