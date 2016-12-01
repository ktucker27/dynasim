package ode;

import utils.DynaComplex;
import coupling.DynaComplexCoupling;

public class CorrelationODEs implements DynaComplexODEs {

    private int n;
    private double gamma;
    private double w;
    private DynaComplexCoupling coupling;
    private double[] d;
    
    DynaComplex[] szs;
    
    DynaComplex[] c1;
    
    DynaComplex t1;
    DynaComplex sum;
    
    /**
     * @param n number of atoms
     * @param gamma spontaneous decay rate
     * @param w pumping rate
     * @param coupling coupling coefficients
     * @param d size n vector of natural frequencies
     * @param szs size n vector of the steady state values of sigma_z
     */
    public CorrelationODEs(int n, double gamma, double w, DynaComplexCoupling coupling, double[] d, DynaComplex[] szs) {
        super();
        this.n = n;
        this.gamma = gamma;
        this.w = w;
        this.coupling = coupling;
        this.d = d;
        
        this.szs = szs;
        
        t1 = new DynaComplex(0, 0);
        
        sum = new DynaComplex(0, 0);
        
        initConstants();
    }
    
    private void initConstants() {
        c1 = new DynaComplex[n];
        for(int i = 0; i < n; ++i) {
            c1[i] = new DynaComplex(-0.5*(gamma+w), -d[i]);
        }
    }
    
    @Override
    public void computeDerivatives(double t, DynaComplex[] z, DynaComplex[] zDot) {
        int idx = 0;
        for(int a = 0; a < n; ++a) {
            for(int b = 0; b < n; ++b) {
//                if(a == b) {
//                    zDot[idx].set(0, 0);
//                    ++idx;
//                    continue;
//                }
                
                sum.set(0,0);
                for(int j = 0; j < n; ++j) {
                    if(j == a || j == b) continue;
                    
                    t1.set(szs[a]).multiply(z[j*n + b]).multiply(coupling.getF(a, j));
                    sum.add(t1);
                }
                sum.multiply(0.5);
                
                zDot[idx].set(c1[a]).multiply(z[idx]);
                zDot[idx].add(t1.set(szs[a]).multiply(z[b*n + b]).multiply(coupling.getF(a, b)).multiply(0.5));
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
