package ode;

import utils.DynaComplex;

public class CorrelationAllToAllODEs implements DynaComplexODEs {

    private int n;
    private double gamma;
    private double w;
    private DynaComplex alpha;
    private double[] d;
    
    DynaComplex[] szs;
    
    DynaComplex[] pmColSums;
    
    DynaComplex[] c1;
    
    DynaComplex t1, t2;
    DynaComplex sum;
    
    /**
     * @param params specifies equation parameters
     * @param szs size n vector of the steady state values of sigma_z
     */
    public CorrelationAllToAllODEs(CumulantParams params, DynaComplex[] szs) {
        this(params.getN(), params.getGamma(), params.getW(), params.getAlpha(), params.getD(), szs);
    }
    
    /**
     * @param n number of atoms
     * @param gamma spontaneous decay rate
     * @param w pumping rate
     * @param coupling coupling coefficients
     * @param d size n vector of natural frequencies
     * @param szs size n vector of the steady state values of sigma_z
     */
    public CorrelationAllToAllODEs(int n, double gamma, double w, DynaComplex alpha, double[] d, DynaComplex[] szs) {
        super();
        this.n = n;
        this.gamma = gamma;
        this.w = w;
        this.alpha = new DynaComplex(alpha);
        this.d = d;
        
        this.szs = szs;
        
        pmColSums = new DynaComplex[n];
        
        t1 = new DynaComplex(0, 0);
        t2 = new DynaComplex(0, 0);
        
        sum = new DynaComplex(0, 0);
        
        initConstants();
    }
    
    private void initConstants() {
        c1 = new DynaComplex[n];
        for(int i = 0; i < n; ++i) {
            c1[i] = new DynaComplex(-0.5*(gamma+w), -d[i]);
            pmColSums[i] = new DynaComplex(0, 0);
        }
    }
    
    @Override
    public void computeDerivatives(double t, DynaComplex[] z, DynaComplex[] zDot) {
        // Precompute the sums.  This will be O(n^2), but by using these in the
        // equations to come, we can avoid O(n^3) operations
        for(int a = 0; a < n; ++a) {
            pmColSums[a].set(0, 0);
            for(int b = 0; b < n; ++b) {
                pmColSums[a].add(z[b*n + a]);
            }
        }
        
        int idx = 0;
        for(int a = 0; a < n; ++a) {
            for(int b = 0; b < n; ++b) {
                sum.set(pmColSums[b]).subtract(z[idx]);
                sum.multiply(szs[a]).multiply(t1.set(alpha).conjugate());
                sum.multiply(gamma*0.5);
                
                zDot[idx].set(c1[a]).multiply(z[idx]);
                zDot[idx].add(sum);
                
                ++idx;
            }
        }
    }

    @Override
    public int getDimension() {
        return n*n;
    }

}
