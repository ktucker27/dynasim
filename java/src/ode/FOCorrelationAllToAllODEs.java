package ode;

import eval.CumulantEval;
import utils.DynaComplex;

public class FOCorrelationAllToAllODEs implements DynaComplexODEs {
    
    /**
     * Fourth order cumulant version of the two-time correlation solver
     */

    private int n;
    private double gamma;
    private double w;
    private DynaComplex alpha;
    private double[] d;
    
    CumulantEval eval;
    DynaComplex[] steady;

    double[] szs;
    double[][] zzs;
    DynaComplex[][] pms;
    
    DynaComplex[][] tripleSums;
    DynaComplex[] ttColSums;
    DynaComplex[] pmRowSums;
    DynaComplex[] pmColSums;
    
    DynaComplex[] c1;
    DynaComplex[] c2;
    
    DynaComplex t1, t2;
    
    /**
     * @param params specifies equation parameters
     * @param z complex output of cumulant expansion equations
     */
    public FOCorrelationAllToAllODEs(CumulantParams params, DynaComplex[] z) {
        this(params.getN(), params.getGamma(), params.getW(), params.getAlpha(), params.getD(), z);
    }
    
    /**
     * @param n number of atoms
     * @param gamma spontaneous decay rate
     * @param w pumping rate
     * @param coupling coupling coefficients
     * @param d size n vector of natural frequencies
     * @param z complex output of cumulant expansion equations
     */
    public FOCorrelationAllToAllODEs(int n, double gamma, double w, DynaComplex alpha, double[] d, DynaComplex[] z) {
        super();
        this.n = n;
        this.gamma = gamma;
        this.w = w;
        this.alpha = new DynaComplex(alpha);
        this.d = d;
        
        eval = new CumulantEval(n);
        this.steady = z;
        
        szs = new double[n];
        zzs = new double[n][n];
        pms = new DynaComplex[n][n];
        
        tripleSums = new DynaComplex[n][n];
        ttColSums = new DynaComplex[n];
        pmRowSums = new DynaComplex[n];
        pmColSums = new DynaComplex[n];
        
        t1 = new DynaComplex(0, 0);
        t2 = new DynaComplex(0, 0);
        
        initConstants();
    }
    
    private void initConstants() {
        c1 = new DynaComplex[n];
        c2 = new DynaComplex[n];

        for(int i = 0; i < n; ++i) {
            ttColSums[i] = new DynaComplex(0, 0);
            pmRowSums[i] = new DynaComplex(0, 0);
            pmColSums[i] = new DynaComplex(0, 0);
            
            for(int j = 0; j < n; ++j) {
                tripleSums[i][j] = new DynaComplex(0, 0);
            }
        }
        
        for(int i = 0; i < n; ++i) {
            c1[i] = new DynaComplex(-0.5*(gamma+w), -d[i]);
            c2[i] = new DynaComplex(-1.5*(gamma+w), -d[i]);
            
            szs[i] = eval.getSingle(2, i, steady, t1).getReal();
            
            for(int j = 0; j < n; ++j) {
                if(i == j) continue;
                
                zzs[i][j] = eval.getDouble(2, 2, i, j, steady, t1).getReal();
                pms[i][j] = new DynaComplex(eval.getDouble(0, 1, i, j, steady, t1));
                
                pmRowSums[i].add(pms[i][j]);
                pmColSums[j].add(pms[i][j]);
            }
        }
    }
    
    private int getTripleIndex(int a, int b, int c) {
        int idx = n*n + c*n*(n-1) + a*(n-1) + b;
        
        if(b > a) --idx;
        
        return idx;
    }
    
    @Override
    public void computeDerivatives(double t, DynaComplex[] z, DynaComplex[] zDot) {
        // Precompute the sums.  This will be O(n^3), but by using these in the
        // equations to come, we can avoid O(n^4) operations
        for(int a = 0; a < n; ++a) {
            ttColSums[a].set(0, 0);
            for(int b = 0; b < n; ++b) {
                ttColSums[a].add(z[b*n + a]);
                
                tripleSums[a][b].set(0, 0);
                for(int j = 0; j < n; ++j) {
                    if(j == a) continue;
                    tripleSums[a][b].add(z[getTripleIndex(a, j, b)]);
                }
            }
        }
        
        int idx = 0;
        
        // <sigma_a^+(t + tau) sigma_b^-(t)>
        for(int a = 0; a < n; ++a) {
            for(int b = 0; b < n; ++b) {
                
                t1.set(tripleSums[a][b]).multiply(t2.set(alpha).conjugate()).multiply(0.5*gamma);
                
                zDot[idx].set(c1[a]).multiply(z[idx]);
                zDot[idx].add(t1);
                
                ++idx;
            }
        }
        
        // <sigma_a^z(t + tau) sigma_b^+(t + tau) sigma_c^-(t)>
        for(int c = 0; c < n; ++c) {
            for(int a = 0; a < n; ++a) {
                for(int b = 0; b < n; ++b) {
                    if(a == b) continue;
                    
                    zDot[idx].set(c2[b]).multiply(z[getTripleIndex(a, b, c)]);
                    zDot[idx].subtract(t1.set(z[b*n+c]).multiply(gamma - w));
                    zDot[idx].subtract(t1.set(z[a*n+c]).multiply(alpha).multiply(0.5*gamma));
                    zDot[idx].subtract(t1.set(z[getTripleIndex(b, a, c)]).multiply(gamma*alpha.getReal()));

                    t1.set(tripleSums[b][c]).subtract(z[getTripleIndex(b, a, c)]).multiply(szs[a]);
                    t1.add(t2.set(tripleSums[a][c]).subtract(z[getTripleIndex(a, b, c)]).multiply(szs[b]));
                    t1.add(t2.set(ttColSums[c]).subtract(z[a*n+c]).subtract(z[b*n+c]).multiply(zzs[a][b] - 2*szs[a]*szs[b]));
                    t1.multiply(t2.set(alpha).conjugate()).multiply(0.5*gamma);
                    zDot[idx].add(t1);

                    t1.set(pmRowSums[b]).subtract(pms[b][a]).multiply(z[a*n+c]);
                    t1.add(t2.set(pmRowSums[a]).subtract(pms[a][b]).multiply(z[b*n+c]));
                    t1.multiply(alpha).multiply(gamma);
                    zDot[idx].subtract(t1);

                    t1.set(ttColSums[c]).subtract(z[a*n+c]).subtract(z[b*n+c]).multiply(pms[b][a]);
                    t1.add(t2.set(pmColSums[a]).subtract(pms[b][a]).multiply(z[b*n+c]));
                    t1.multiply(t2.set(alpha).conjugate()).multiply(gamma);
                    zDot[idx].subtract(t1);
                    
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
