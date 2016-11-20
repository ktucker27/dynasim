package ode;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;

import utils.DynaComplex;
import coupling.DynaComplexCoupling;

public class DynaCumulantODEs implements DynaComplexODEs {
    
    private int n;
    private double gamma;
    private double w;
    private DynaComplexCoupling coupling;
    private double[] d;

    double diff_gw;
    DynaComplex migamma;
    private DynaComplex[] c1;
    private DynaComplex[] c2;
    private DynaComplex[][] c3;
    private DynaComplex[][] c4;
    
    private DynaComplex t1, t2;
    private DynaComplex ct1, ct2, ct3;
    private DynaComplex sum1, sum2, sum3, sum4;
    
    int[] startIdx;
    
    /**
     * @param n number of oscillators
     * @param gamma spontaneous decay rate
     * @param w pumping rate
     * @param feff effective mean field coupling
     * @param d size n vector of natural frequencies
     */
    public DynaCumulantODEs(int n, double gamma, double w, DynaComplexCoupling coupling, double[] d) {
        super();
        this.n = n;
        this.gamma = gamma;
        this.w = w;
        this.coupling = coupling;
        this.d = d;
        
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
        sum4 = new DynaComplex(0, 0);
        
        t1 = new DynaComplex(0, 0);
        t2 = new DynaComplex(0, 0);

        ct1 = new DynaComplex(0, 0);
        ct2 = new DynaComplex(0, 0);
        ct3 = new DynaComplex(0, 0);

        initConstants();
    }
    
    private void initConstants() {
        diff_gw = gamma - w;
        migamma = new DynaComplex(0,-gamma);
        c1 = new DynaComplex[n];
        for(int i = 0; i < n; ++i) {
            c1[i] = new DynaComplex(-0.5*(gamma + w), -d[i]);
        }
        
        c2 = new DynaComplex[n];
        for(int i = 0; i < n; ++i) {
            c2[i] = new DynaComplex(-1.5*(gamma + w), -d[i]);
        }
        
        c3 = new DynaComplex[n][n];
        for(int i = 0; i < n; ++i) {
            for(int j = i+1; j < n; ++j) {
                if(i == j) continue;
                
                c3[i][j] = new DynaComplex(-(gamma + w), -(d[i] - d[j]));
            }
        }
        
        c4 = new DynaComplex[n][n];
        for(int i = 0; i < n; ++i) {
            for(int j = i+1; j < n; ++j) {
                if(i == j) continue;
                
                c4[i][j] = new DynaComplex(-(gamma + w), -(d[i] + d[j]));
            }
        }
    }
    
    private int getTriIdx(int i, int j) {
        if(i > j) {
            int t = j;
            j = i;
            i = t;
        }
        
        //return n*(n-1)/2 - (n-i)*(n-i-1)/2 + j-1;
        return n*i - i*(i+1)/2 + j - i - 1; // Simplified
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
        switch(al) {
        case 0:
            return ans.set(z[a]);
        case 1:
            return ans.set(z[a]).conjugate();
        case 2:
            return ans.set(z[startIdx[2] + a]);
        }
        
        // This should never happen
        System.err.println("Invalid superscript passed to cumulantSingle: " + al);
        
        return null;
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
    
    private DynaComplex cumulant(int al, int bt, int gm, int a, int b, int c, DynaComplex[] z, DynaComplex ans) {
        cumulantDouble(al, bt, a, b, z, ans).multiply(cumulantSingle(gm, c, z, ct1));
        ans.add(cumulantDouble(bt, gm, b, c, z, ct1).multiply(cumulantSingle(al, a, z, ct2)));
        ans.add(cumulantDouble(al, gm, a, c, z, ct1).multiply(cumulantSingle(bt, b, z, ct2)));
        ans.subtract(cumulantSingle(al, a, z, ct1).multiply(cumulantSingle(bt, b, z, ct2)).multiply(cumulantSingle(gm, c, z, ct3)).multiply(2.0));
        
        return ans;
    }

    @Override
    public void computeDerivatives(double t, DynaComplex[] z, DynaComplex[] zDot)
            throws MaxCountExceededException, DimensionMismatchException {
        int idx;
        
        // sigma_a^+
        for(int a = 0; a < n; ++a) {
            sum1.set(0, 0);
            idx = startIdx[1] + a*(n-1);
            for(int b = 0; b < n; ++b) {
                if(a == b) continue;
                
                sum1.add(t1.set(coupling.getAlpha(a, b)).conjugate().multiply(z[idx]));
                ++idx;
            }
            sum1.multiply(0.5*gamma);
            
            zDot[a].set(z[a]).multiply(c1[a]).add(sum1);
        }
        
        // sigma_a^z sigma_b^+
        idx = startIdx[1];
        for(int a = 0; a < n; ++a) {
            for(int b = 0; b < n; ++b) {
                if(a == b) continue;
                
                sum1.set(0, 0);
                for(int j = 0; j < n; ++j) {
                    if(j == a || j == b) continue;
                    
                    sum1.add(t1.set(coupling.getAlpha(b, j)).conjugate().multiply(cumulant(2, 2, 0, a, b, j, z, t2)));
                }
                sum1.multiply(0.5*gamma);
                
                sum2.set(0, 0);
                for(int j = 0; j < n; ++j) {
                    if(j == a || j == b) continue;
                    
                    sum2.add(t1.set(coupling.getAlpha(a, j)).multiply(cumulant(0, 0, 1, a, b, j, z, t2)));
                }
                sum2.multiply(-1.0*gamma);
                
                sum3.set(0, 0);
                for(int j = 0; j < n; ++j) {
                    if(j == a || j == b) continue;
                    
                    sum3.add(t1.set(coupling.getAlpha(a, j)).conjugate().multiply(cumulant(1, 0, 0, a, b, j, z, t2)));
                }
                sum3.multiply(-1.0*gamma);
                
                zDot[idx].set(z[idx]).multiply(c2[b]).add(t1.set(z[b]).multiply(-diff_gw)).add(t1.set(coupling.getAlpha(a, b)).multiply(z[a]).multiply(-0.5*gamma));
                zDot[idx].add(t1.set(z[n+getRecIdx(b,a)]).multiply(coupling.getF(a, b)).multiply(-gamma));
                zDot[idx].add(sum1).add(sum2).add(sum3);
                ++idx;
            }
        }
        
        // sigma_a^z
        for(int a = 0; a < n; ++a) {
            sum1.set(0, 0);
            for(int b = 0; b < n; ++b) {
                if(a == b) continue;
                
                sum1.add(cumulantDouble(0, 1, a, b, z, t1).subtract(cumulantDouble(1, 0, a, b, z, t2)).multiply(coupling.getG(a, b)));
            }
            sum1.multiply(migamma);
            
            sum2.set(0, 0);
            for(int b = 0; b < n; ++b) {
                if(a == b) continue;
                
                sum2.add(cumulantDouble(0, 1, a, b, z, t1).add(cumulantDouble(1, 0, a, b, z, t2)).multiply(coupling.getF(a, b)));
            }
            sum2.multiply(-gamma);
            
            zDot[startIdx[2] + a].set(z[startIdx[2] + a]).multiply(-(gamma + w)).add(w - gamma).add(sum1).add(sum2);
        }
        
        // sigma_a^+ sigma_b^-
        idx = startIdx[3];
        for(int a = 0; a < n; ++a) {
            for(int b = a+1; b < n; ++b) {
                sum1.set(0, 0);
                for(int j = 0; j < n; ++j) {
                    if(j == a || j == b) continue;
                    
                    sum1.add(t1.set(coupling.getAlpha(a, j)).conjugate().multiply(cumulant(2, 1, 0, a, b, j, z, t2)));
                }
                sum1.multiply(0.5*gamma);
                
                sum2.set(0, 0);
                for(int j = 0; j < n; ++j) {
                    if(j == a || j == b) continue;
                    
                    sum2.add(t1.set(coupling.getAlpha(b, j)).multiply(cumulant(2, 0, 1, b, a, j, z, t2)));
                }
                sum2.multiply(0.5*gamma);
                
                zDot[idx].set(z[idx]).multiply(c3[a][b]).add(t1.set(migamma).multiply(0.25).multiply(coupling.getG(a,b)).multiply(t2.set(z[startIdx[2] + a]).subtract(z[startIdx[2] + b])));
                zDot[idx].add(cumulantDouble(2, 2, a, b, z, t1).add((t2.set(z[startIdx[2] + a]).add(z[startIdx[2] + b])).multiply(0.5)).multiply(0.5*gamma*coupling.getF(a,b)));
                zDot[idx].add(sum1).add(sum2);
                ++idx;
            }
        }
        
        // sigma_a^z sigma_b^z
        idx = startIdx[4];
        for(int a = 0; a < n; ++a) {
            for(int b = a+1; b < n; ++b) {
                sum1.set(0, 0);
                for(int j = 0; j < n; ++j) {
                    if(j == a || j == b) continue;
                    
                    sum1.add(cumulant(0, 2, 1, a, b, j, z, t1).subtract(cumulant(1, 2, 0, a, b, j, z, t2)).multiply(coupling.getG(a, j)));
                }
                sum1.multiply(migamma);
                
                sum2.set(0, 0);
                for(int j = 0; j < n; ++j) {
                    if(j == a || j == b) continue;
                    
                    sum2.add(cumulant(0, 2, 1, b, a, j, z, t1).subtract(cumulant(1, 2, 0, b, a, j, z, t2)).multiply(coupling.getG(b, j)));
                }
                sum2.multiply(migamma);
                
                sum3.set(0, 0);
                for(int j = 0; j < n; ++j) {
                    if(j == a || j == b) continue;
                    
                    sum3.add(cumulant(0, 2, 1, a, b, j, z, t1).add(cumulant(1, 2, 0, a, b, j, z, t2)).multiply(coupling.getF(a, j)));
                }
                sum3.multiply(-gamma);
                
                sum4.set(0, 0);
                for(int j = 0; j < n; ++j) {
                    if(j == a || j == b) continue;
                    
                    sum4.add(cumulant(0, 2, 1, b, a, j, z, t1).add(cumulant(1, 2, 0, b, a, j, z, t2)).multiply(coupling.getF(b, j)));
                }
                sum4.multiply(-gamma);
                
                zDot[idx].set(z[startIdx[2] + a]).add(z[startIdx[2] + b]).multiply(w - gamma);
                zDot[idx].subtract(cumulantDouble(2, 2, a, b, z, t1).multiply(2*(gamma + w)));
                zDot[idx].add(cumulantDouble(0, 1, a, b, z, t1).add(cumulantDouble(1, 0, a, b, z, t2)).multiply(2*gamma*coupling.getF(a, b)));
                zDot[idx].add(sum1).add(sum2).add(sum3).add(sum4);
                ++idx;
            }
        }
        
        // sigma_a^+ sigma_b^+
        idx = startIdx[5];
        for(int a = 0; a < n; ++a) {
            for(int b = a+1; b < n; ++b) {
                sum1.set(0, 0);
                for(int j = 0; j < n; ++j) {
                    if(j == a || j == b) continue;
                    
                    sum1.add(cumulant(2, 0, 0, a, b, j, z, t1).multiply(t2.set(coupling.getAlpha(a, j)).conjugate()));
                }
                sum1.multiply(0.5*gamma);
                
                sum2.set(0, 0);
                for(int j = 0; j < n; ++j) {
                    if(j == a || j == b) continue;
                    
                    sum2.add(cumulant(2, 0, 0, b, a, j, z, t1).multiply(t2.set(coupling.getAlpha(b, j)).conjugate()));
                }
                sum2.multiply(0.5*gamma);
                
                zDot[idx].set(z[idx]).multiply(c4[a][b]).add(sum1).add(sum2);
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
