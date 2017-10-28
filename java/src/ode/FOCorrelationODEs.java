package ode;

import utils.DynaComplex;
import coupling.DynaComplexCoupling;

public class FOCorrelationODEs implements DynaComplexODEs {
    
    /**
     * Fourth order cumulant version of the two-time correlation solver
     */

    private int n;
    private double gamma;
    private double w;
    private DynaComplexCoupling coupling;
    private double[] d;
    
    DynaComplex[] szs;
    
    DynaComplex[] c1;
    DynaComplex[] c2;
    
    DynaComplex t1, t2;
    DynaComplex sum;
    
    /**
     * @param n number of atoms
     * @param gamma spontaneous decay rate
     * @param w pumping rate
     * @param coupling coupling coefficients
     * @param d size n vector of natural frequencies
     * @param szs size n vector of the steady state values of sigma_z
     */
    public FOCorrelationODEs(int n, double gamma, double w, DynaComplexCoupling coupling, double[] d, DynaComplex[] szs) {
        super();
        this.n = n;
        this.gamma = gamma;
        this.w = w;
        this.coupling = coupling;
        this.d = d;
        
        this.szs = szs;
        
        t1 = new DynaComplex(0, 0);
        t2 = new DynaComplex(0, 0);
        
        sum = new DynaComplex(0, 0);
        
        initConstants();
    }
    
    private void initConstants() {
        c1 = new DynaComplex[n];
        c2 = new DynaComplex[n];
        for(int i = 0; i < n; ++i) {
            c1[i] = new DynaComplex(-0.5*(gamma+w), -d[i]);
            c2[i] = new DynaComplex(-1.5*(gamma+w), -d[i]);
        }
    }
    
    private int getTripleIndex(int a, int b, int c) {
        int idx = n*n + c*n*(n-1) + a*(n-1) + b;
        
        if(b > a) --idx;
        
        return idx;
    }
    
    @Override
    public void computeDerivatives(double t, DynaComplex[] z, DynaComplex[] zDot) {
        int idx = 0;
        
        // <sigma_a^+(t + tau) sigma_b^-(t)>
        for(int a = 0; a < n; ++a) {
            for(int b = 0; b < n; ++b) {
                
                sum.set(0,0);
                for(int j = 0; j < n; ++j) {
                    if(j == a) continue;
                    
                    t1.set(z[getTripleIndex(a,j,b)]).multiply(t2.set(coupling.getAlpha(a, j)).conjugate());
                    sum.add(t1);
                }
                sum.multiply(gamma*0.5);
                
                zDot[idx].set(c1[a]).multiply(z[idx]);
                zDot[idx].add(sum);
                
                ++idx;
            }
        }
        
        // <sigma_a^z(t + tau) sigma_b^+(t + tau) sigma_c^-(t)>
        for(int c = 0; c < n; ++c) {
            for(int a = 0; a < n; ++a) {
                for(int b = 0; b < n; ++b) {
                    if(a == b) continue;
                    
                    sum.set(0,0);
                    for(int j = 0; j < n; ++j) {
                        if(j == a || j == b) continue;
                        
                        t1.set(z[getTripleIndex(b, j, c)]).multiply(szs[a]);
                        t1.add(t2.set(z[getTripleIndex(a, j, c)]).multiply(szs[b]));
                        t1.subtract(t2.set(z[j*n+c]).multiply(szs[a]).multiply(szs[b]));
                        t1.multiply(t2.set(coupling.getAlpha(b, j)).conjugate());
                        sum.add(t1);
                    }
                    sum.multiply(0.5*gamma);
                    
                    zDot[idx].set(c2[b]).multiply(z[getTripleIndex(a, b, c)]);
                    zDot[idx].subtract(t1.set(z[b*n+c]).multiply(gamma - w));
                    zDot[idx].subtract(t1.set(z[a*n+c]).multiply(coupling.getAlpha(a, b)).multiply(0.5*gamma));
                    zDot[idx].subtract(t1.set(z[getTripleIndex(b, a, c)]).multiply(gamma*coupling.getAlpha(a, b).getReal()));
                    zDot[idx].add(sum);
                    
                    ++idx;
                }
            }
        }
    }

    @Override
    public int getDimension() {
        return n*n + n*n*(n-1);
    }

}
