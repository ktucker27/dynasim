package coupling;

import org.apache.commons.math3.complex.Complex;

public class RingCoupling implements ComplexCoupling {

    private double[] f;
    private double[] g;
    private Complex[] alpha;
    
    /**
     * @param n the number of oscillators
     * @param a angular separation
     */
    public RingCoupling(int n, double a) {
        super();
        
        f = new double[n];
        g = new double[n];
        alpha = new Complex[n];
        
        f[0] = 0.0;
        g[0] = 0.0;
        alpha[0] = new Complex(0.0);

        for(int m = 1; m < n; ++m) {
            double am = a*m;
            double am2 = am*am;
            double am3 = am2*am;
            
            double sinAm = Math.sin(am);
            double cosAm = Math.cos(am);
            
            f[m] = 1.5*(sinAm/am + cosAm/am2 - sinAm/am3);
            g[m] = 1.5*(-cosAm/am + sinAm/am2 + cosAm/am3);
            //g[m] = 0.0;
            alpha[m] = new Complex(f[m], g[m]);
        }
    }

    @Override
    public double getF(int row, int col) {
        int m = Math.abs(col - row);
        return f[m];
    }

    @Override
    public double getG(int row, int col) {
        int m = Math.abs(col - row);
        return g[m];
    }
    
    @Override
    public Complex getAlpha(int row, int col) {
        int m = Math.abs(col - row);
        return alpha[m];
    }
}
