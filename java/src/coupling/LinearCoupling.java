package coupling;

import utils.DynaComplex;

public class LinearCoupling implements DynaComplexCoupling {

    double[] fs;
    DynaComplex[] alphas;
    
    public LinearCoupling(double al, double a, int n) {
        fs = new double[n];
        alphas = new DynaComplex[n];
        fs[0] = 0.0;
        alphas[0] = new DynaComplex(0, 0);
        for(int i = 1; i < n; ++i) {
            fs[i] = Math.pow(Math.abs(i)*a, -al);
            alphas[i] = new DynaComplex(fs[i], 0.0);
            //System.out.println(fs[i]);
        }
    }
    
    @Override
    public double getF(int row, int col) {
        return fs[Math.max(row - col, col - row)];
    }

    @Override
    public double getG(int row, int col) {
        return 0;
    }

    @Override
    public DynaComplex getAlpha(int row, int col) {
        return alphas[Math.max(row - col, col - row)];
    }

}
