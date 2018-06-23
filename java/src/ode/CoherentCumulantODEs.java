package ode;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;

import utils.DynaComplex;

public class CoherentCumulantODEs implements DynaComplexODEs {
    
    private int n;
    private double gamma;
    private double w;
    private double gel;
    private DynaComplex alpha;
    private DynaComplex alpha_bar;
    private double[] d;
    private double omega;

    double diff_gw;
    DynaComplex migamma;
    DynaComplex miomd2;
    private DynaComplex[] c1;
    private DynaComplex[] c2;
    private DynaComplex[][] c3;
    private DynaComplex[][] c4;
    
    private DynaComplex t1, t2;
    private DynaComplex ct1, ct2, ct3;
    private DynaComplex sum1, sum2, sum3;
    
    private DynaComplex psum, zsum;
    private DynaComplex[] zp_row_sum, zp_col_sum, pmsum, zzsum, ppsum;
    
    int[] startIdx;
    
    /**
     * @param params input object specifying parameters
     */
    public CoherentCumulantODEs(SystemParams params, double omega) {
        super();
        this.n = params.getN();
        this.gamma = params.getGamma();
        this.w = params.getW();
        this.gel = params.getGel();
        this.alpha = new DynaComplex(params.getAlpha());
        this.alpha_bar = new DynaComplex(alpha.getReal(), -alpha.getImaginary());
        this.d = params.getD();
        this.omega = omega;
        
        init();
    }
    
    /**
     * @param n number of atoms
     * @param gamma spontaneous decay rate
     * @param w pumping rate
     * @param feff effective mean field coupling
     * @param d size n vector of natural frequencies
     */
    public CoherentCumulantODEs(int n, double gamma, double w, double gel, DynaComplex alpha, double[] d, double omega) {
        super();
        this.n = n;
        this.gamma = gamma;
        this.w = w;
        this.gel = gel;
        this.alpha = new DynaComplex(alpha.getReal(), alpha.getImaginary());
        this.alpha_bar = new DynaComplex(alpha.getReal(), -alpha.getImaginary());
        this.d = d;
        this.omega = omega;
        
        init();
    }
    
    private void init() {
        startIdx = new int[6];
        startIdx[0] = 0;
        startIdx[1] = n;
        startIdx[2] = startIdx[1] + n*(n-1);
        startIdx[3] = startIdx[2] + n;
        startIdx[4] = startIdx[3] + n*(n-1)/2;
        startIdx[5] = startIdx[4] + n*(n-1)/2;
        
        sum1 = new DynaComplex(0, 0);
        sum2 = new DynaComplex(0, 0);
        sum3 = new DynaComplex(0, 0);
        
        t1 = new DynaComplex(0, 0);
        t2 = new DynaComplex(0, 0);

        ct1 = new DynaComplex(0, 0);
        ct2 = new DynaComplex(0, 0);
        ct3 = new DynaComplex(0, 0);
        
        psum = new DynaComplex(0, 0);
        zsum = new DynaComplex(0, 0);
        
        zp_row_sum = new DynaComplex[n];
        zp_col_sum = new DynaComplex[n];
        pmsum = new DynaComplex[n];
        zzsum = new DynaComplex[n];
        ppsum = new DynaComplex[n];
        for(int i = 0; i < n; ++i) {
            zp_row_sum[i] = new DynaComplex(0, 0);
            zp_col_sum[i] = new DynaComplex(0, 0);
            pmsum[i] = new DynaComplex(0, 0);
            zzsum[i] = new DynaComplex(0, 0);
            ppsum[i] = new DynaComplex(0, 0);
        }
        
        initConstants();
    }
    
    private void initConstants() {
        diff_gw = gamma - w;
        migamma = new DynaComplex(0,-gamma);
        miomd2 = new DynaComplex(0,-0.5*omega);
        c1 = new DynaComplex[n];
        for(int i = 0; i < n; ++i) {
            c1[i] = new DynaComplex(-0.5*(gamma + w + gel), -d[i]);
        }
        
        c2 = new DynaComplex[n];
        for(int i = 0; i < n; ++i) {
            c2[i] = new DynaComplex(-1.5*(gamma + w) - 0.5*gel, -d[i]);
        }
        
        c3 = new DynaComplex[n][n];
        for(int i = 0; i < n; ++i) {
            for(int j = i+1; j < n; ++j) {
                if(i == j) continue;
                
                c3[i][j] = new DynaComplex(-(gamma + w + gel), -(d[i] - d[j]));
            }
        }
        
        c4 = new DynaComplex[n][n];
        for(int i = 0; i < n; ++i) {
            for(int j = i+1; j < n; ++j) {
                if(i == j) continue;
                
                c4[i][j] = new DynaComplex(-(gamma + w + gel), -(d[i] + d[j]));
            }
        }
    }
    
    public void setW(double w) {
        this.w = w;
        initConstants();
    }
    
    private int getTriIdx(int i, int j) {
        if(i > j) {
            int t = j;
            j = i;
            i = t;
        }
        
        return n*i - i*(i+1)/2 + j - i - 1;
    }
    
    private int getRecIdx(int i, int j) {
        if(j < i) {
            return i*(n-1) + j;
        }
        
        return i*(n-1) + j - 1;
    }
    
    // In the following, for the superscript indices:
    // 0 = +, 1 = -, 2 = z
    private DynaComplex cumulantSingle(int al, int a, DynaComplex[] z, DynaComplex ans) {
        ans.set(z[a + startIdx[2]*(al/2)]);
        if(al == 1) ans.conjugate();
        return ans;
    }
    
    private DynaComplex cumulantSingleSum(int al, DynaComplex ans) {
        switch(al) {
        case 0:
            ans.set(psum);
            break;
        case 1:
            ans.set(psum).conjugate();
            break;
        case 2:
            ans.set(zsum);
        }
        
        return ans;
    }

    private DynaComplex cumulantDouble(int al, int bt, int a, int b, DynaComplex[] z, DynaComplex ans) {
        switch(al) {
        case 0:
            switch(bt) {
            case 0:
                return ans.set(z[startIdx[5] + getTriIdx(a,b)]);
            case 1:
                ans.set(z[startIdx[3] + getTriIdx(a,b)]);
                if(a > b) return ans.conjugate();
                return ans;
            case 2:
                return ans.set(z[startIdx[1] + getRecIdx(b,a)]);
            }
        case 1:
            switch(bt) {
            case 0:
                ans.set(z[startIdx[3] + getTriIdx(a,b)]);
                if(a < b) return ans.conjugate();
                return ans;
            case 1:
                return ans.set(z[startIdx[5] + getTriIdx(a,b)]).conjugate();
            case 2:
                return ans.set(z[startIdx[1] + getRecIdx(b,a)]).conjugate();
            }
        case 2:
            switch(bt) {
            case 0:
                return ans.set(z[startIdx[1] + getRecIdx(a,b)]);
            case 1:
                return ans.set(z[startIdx[1] + getRecIdx(a,b)]).conjugate();
            case 2:
                return ans.set(z[startIdx[4] + getTriIdx(a,b)]);
            }
        }
        
        // This should never happen
        System.err.println("Invalid superscripts passed to cumulantDouble: " + al + " " + bt);
        
        return null;
    }
    
    private DynaComplex cumulantDoubleSum(int al, int bt, int a, DynaComplex ans) {
        switch(al) {
        case 0:
            switch(bt) {
            case 0:
                ans.set(ppsum[a]);
                break;
            case 1:
                ans.set(pmsum[a]);
                break;
            case 2:
                ans.set(zp_col_sum[a]);
                break;
            }
            break;
        case 1:
            switch(bt) {
            case 0:
                ans.set(pmsum[a]).conjugate();
                break;
            case 1:
                ans.set(ppsum[a]).conjugate();
                break;
            case 2:
                ans.set(zp_col_sum[a]).conjugate();
                break;
            }
            break;
        case 2:
            switch(bt) {
            case 0:
                ans.set(zp_row_sum[a]);
                break;
            case 1:
                ans.set(zp_row_sum[a]).conjugate();
                break;
            case 2:
                ans.set(zzsum[a]);
                break;
            }
            break;
        }
        return ans;
    }
    
    private DynaComplex cumulantSum(int al, int bt, int gm, int a, int b, DynaComplex[] z, DynaComplex ans) {
        cumulantSingleSum(gm, ans).subtract(cumulantSingle(gm, a, z, ct1)).subtract(cumulantSingle(gm, b, z, ct2));
        ans.multiply(cumulantDouble(al, bt, a, b, z, ct1).subtract(cumulantSingle(al, a, z, ct2).multiply(cumulantSingle(bt, b, z, ct3)).multiply(2.0)));
        ans.add(cumulantSingle(al, a, z, ct1).multiply(cumulantDoubleSum(bt, gm, b, ct2).subtract(cumulantDouble(bt, gm, b, a, z, ct3))));
        ans.add(cumulantSingle(bt, b, z, ct1).multiply(cumulantDoubleSum(al, gm, a, ct2).subtract(cumulantDouble(al, gm, a, b, z, ct3))));
        return ans;
    }

    @Override
    public void computeDerivatives(double t, DynaComplex[] z, DynaComplex[] zDot)
            throws MaxCountExceededException, DimensionMismatchException {
        int idx;
        
        // Start by computing the sums that we need.  This is O(n^2)
        psum.set(0, 0);
        zsum.set(0, 0);
        for(int a = 0; a < n; ++a) {
            psum.add(z[a]);
            zsum.add(z[startIdx[2] + a]);
            
            zp_row_sum[a].set(0, 0);
            zp_col_sum[a].set(0, 0);
            pmsum[a].set(0, 0);
            zzsum[a].set(0, 0);
            ppsum[a].set(0, 0);
            for(int b = 0; b < n; ++b) {
                if(a == b) continue;
                
                zp_row_sum[a].add(z[startIdx[1] + getRecIdx(a, b)]);
                zp_col_sum[a].add(z[startIdx[1] + getRecIdx(b, a)]);
                idx = getTriIdx(a,b);
                zzsum[a].add(z[startIdx[4] + idx]);
                ppsum[a].add(z[startIdx[5] + idx]);
                if(b > a) {
                    pmsum[a].add(z[startIdx[3] + idx]);
                } else {
                    pmsum[a].add(t1.set(z[startIdx[3] + idx]).conjugate());
                }
            }
        }
        
        // sigma_a^+
        for(int a = 0; a < n; ++a) {
            sum1.set(0, 0);
            idx = startIdx[1] + a*(n-1);
            for(int b = 0; b < n; ++b) {
                if(a == b) continue;
                
                sum1.add(t1.set(alpha_bar).multiply(z[idx]));
                ++idx;
            }
            sum1.multiply(0.5*gamma);
            
            zDot[a].set(z[a]).multiply(c1[a]).add(sum1);
            
            // Coherent pumping
            zDot[a].add(t1.set(z[startIdx[2] + a]).multiply(miomd2));
        }
        
        // sigma_a^z sigma_b^+
        idx = startIdx[1];
        for(int a = 0; a < n; ++a) {
            for(int b = 0; b < n; ++b) {
                if(a == b) continue;
                
                // z z +
                cumulantSum(2, 2, 0, a, b, z, sum1);
                sum1.multiply(alpha_bar);
                sum1.multiply(0.5*gamma);
                
                // + + -
                cumulantSum(0, 0, 1, a, b, z, sum2);
                sum2.multiply(alpha);
                sum2.multiply(-gamma);
                
                // - + +
                cumulantSum(1, 0, 0, a, b, z, sum3);
                sum3.multiply(alpha_bar);
                sum3.multiply(-gamma);
                
                zDot[idx].set(z[idx]).multiply(c2[b]).add(t1.set(z[b]).multiply(-diff_gw)).add(t1.set(alpha).multiply(z[a]).multiply(-0.5*gamma));
                zDot[idx].add(t1.set(z[n+getRecIdx(b,a)]).multiply(alpha.getReal()).multiply(-gamma));
                zDot[idx].add(sum1).add(sum2).add(sum3);
                
                // Coherent pumping
                cumulantDouble(0, 0, a, b, z, t1);
                cumulantDouble(1, 0, a, b, z, t2);
                
                zDot[idx].add(t1.subtract(t2).multiply(miomd2).multiply(2));

                cumulantDouble(2, 2, a, b, z, t1);
                
                zDot[idx].add(t1.multiply(miomd2));

                ++idx;
            }
        }
        
        // sigma_a^z
        for(int a = 0; a < n; ++a) {
            sum1.set(0, 0);
            for(int b = 0; b < n; ++b) {
                if(a == b) continue;
                
                t1.set(z[startIdx[3] + getTriIdx(a,b)]);
                if(a > b) {
                    t1.conjugate();
                }
                t1.multiply(alpha).multiply(2.0);
                t1.setImaginary(0);
                sum1.add(t1);
            }
            sum1.multiply(-1.0*gamma);
            
            zDot[startIdx[2] + a].set(z[startIdx[2] + a]).multiply(-(gamma + w)).add(w - gamma).add(sum1);

            // Coherent pumping
            zDot[startIdx[2] + a].add(t1.set(z[a]).subtract(t2.set(z[a]).conjugate()).multiply(miomd2).multiply(2));
        }
        
        // sigma_a^+ sigma_b^-
        idx = startIdx[3];
        for(int a = 0; a < n; ++a) {
            for(int b = a+1; b < n; ++b) {
                // z - +
                cumulantSum(2, 1, 0, a, b, z, sum1);
                sum1.multiply(alpha_bar);
                sum1.multiply(0.5*gamma);
                
                // z + -; b a j
                cumulantSum(2, 0, 1, b, a, z, sum2);
                sum2.multiply(alpha);
                sum2.multiply(0.5*gamma);
                
                zDot[idx].set(z[idx]).multiply(c3[a][b]).add(t1.set(migamma).multiply(0.25).multiply(alpha.getImaginary()).multiply(t2.set(z[startIdx[2] + a]).subtract(z[startIdx[2] + b])));
                zDot[idx].add(t1.set(z[startIdx[4] + getTriIdx(a,b)]).add((t2.set(z[startIdx[2] + a]).add(z[startIdx[2] + b])).multiply(0.5)).multiply(0.5*gamma*alpha.getReal()));
                zDot[idx].add(sum1).add(sum2);
                
                // Coherent pumping
                cumulantDouble(2, 1, a, b, z, t1);
                cumulantDouble(0, 2, a, b, z, t2);
                
                zDot[idx].add(t1.subtract(t2).multiply(miomd2));
                
                ++idx;
            }
        }
        
        // sigma_a^z sigma_b^z
        idx = startIdx[4];
        for(int a = 0; a < n; ++a) {
            for(int b = a+1; b < n; ++b) {
                // + z -
                cumulantSum(0, 2, 1, a, b, z, sum1);
                sum1.multiply(alpha);
                sum1.multiply(-2.0*gamma);
                sum1.setImaginary(0.0);
                
                // + z -; b a j
                cumulantSum(0, 2, 1, b, a, z, sum2);
                sum2.multiply(alpha);
                sum2.multiply(-2.0*gamma);
                sum2.setImaginary(0.0);
                
                int triab = getTriIdx(a,b);
                zDot[idx].set(z[startIdx[3] + triab]).multiply(4.0*gamma*alpha.getReal());
                zDot[idx].setImaginary(0);
                zDot[idx].add(t1.set(z[startIdx[2] + a]).add(z[startIdx[2] + b]).multiply(w - gamma));
                zDot[idx].subtract(t1.set(z[startIdx[4] + triab]).multiply(2*(gamma + w)));
                zDot[idx].add(sum1).add(sum2);
                
                // Coherent pumping
                cumulantDouble(0, 2, a, b, z, t1);
                cumulantDouble(1, 2, a, b, z, t2);
                t1.subtract(t2);
                cumulantDouble(2, 0, a, b, z, t2);
                t1.add(t2);
                cumulantDouble(2, 1, a, b, z, t2);
                t1.subtract(t2).multiply(miomd2).multiply(2);
                
                zDot[idx].add(t1);
                
                ++idx;
            }
        }
        
        // sigma_a^+ sigma_b^+
        idx = startIdx[5];
        for(int a = 0; a < n; ++a) {
            for(int b = a+1; b < n; ++b) {
                // z + +
                cumulantSum(2, 0, 0, a, b, z, sum1);
                sum1.multiply(alpha_bar);
                sum1.multiply(0.5*gamma);
                
                // z + +; b a j
                cumulantSum(2, 0, 0, b, a, z, sum2);
                sum2.multiply(alpha_bar);
                sum2.multiply(0.5*gamma);
                
                zDot[idx].set(z[idx]).multiply(c4[a][b]).add(sum1).add(sum2);
                
                // Coherent pumping
                cumulantDouble(2, 0, a, b, z, t1);
                cumulantDouble(0, 2, a, b, z, t2);
                
                zDot[idx].add(t1.add(t2).multiply(miomd2));
                
                ++idx;
            }
        }
    }

    @Override
    public int getDimension() {
        return n*(n-1) + 3*n*(n-1)/2 + 2*n;
    }
    
    public int getStartIdx(int idx) {
        return 2*startIdx[idx];
    }
}
