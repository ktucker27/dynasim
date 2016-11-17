package ode;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;

import coupling.ComplexCoupling;

public class SynchCumulantlODEs implements ComplexODEs {
    
    private int n;
    private double gamma;
    private double w;
    private ComplexCoupling coupling;
    private double[] d;

    double diff_gw;
    Complex migamma;
    private Complex[] c1;
    private Complex[] c2;
    private Complex[][] c3;
    private Complex[][] c4;
    
    int[] startIdx;
    
    /**
     * @param n number of oscillators
     * @param gamma spontaneous decay rate
     * @param w pumping rate
     * @param feff effective mean field coupling
     * @param d size n vector of natural frequencies
     */
    public SynchCumulantlODEs(int n, double gamma, double w, ComplexCoupling coupling, double[] d) {
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
        
        initConstants();
    }
    
    private void initConstants() {
        diff_gw = gamma - w;
        migamma = new Complex(0,-gamma);
        c1 = new Complex[n];
        for(int i = 0; i < n; ++i) {
            c1[i] = new Complex(-0.5*(gamma + w), -d[i]);
        }
        
        c2 = new Complex[n];
        for(int i = 0; i < n; ++i) {
            c2[i] = new Complex(-1.5*(gamma + w), -d[i]);
        }
        
        c3 = new Complex[n][n];
        for(int i = 0; i < n; ++i) {
            for(int j = i+1; j < n; ++j) {
                if(i == j) continue;
                
                c3[i][j] = new Complex(-(gamma + w), -(d[i] - d[j]));
            }
        }
        
        c4 = new Complex[n][n];
        for(int i = 0; i < n; ++i) {
            for(int j = i+1; j < n; ++j) {
                if(i == j) continue;
                
                c4[i][j] = new Complex(-(gamma + w), -(d[i] + d[j]));
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
    private Complex cumulantSingle(int al, int a, Complex[] z) {
        switch(al) {
        case 0:
            return z[a];
        case 1:
            return z[a].conjugate();
        case 2:
            return z[startIdx[2] + a];
        }
        
        // This should never happen
        System.err.println("Invalid superscript passed to cumulantSingle: " + al);
        
        return null;
    }

    private Complex cumulantDouble(int al, int bt, int a, int b, Complex[] z) {
        Complex ans;
        switch(al) {
        case 0:
            switch(bt) {
            case 0:
                return z[startIdx[5] + getTriIdx(a,b)];
            case 1:
                ans = z[startIdx[3] + getTriIdx(a,b)];
                if(a > b) return ans.conjugate();
                return ans;
            case 2:
                return z[startIdx[1] + getRecIdx(b,a)];
            }
        case 1:
            switch(bt) {
            case 0:
                ans = z[startIdx[3] + getTriIdx(a,b)];
                if(a < b) return ans.conjugate();
                return ans;
            case 1:
                return z[startIdx[5] + getTriIdx(a,b)].conjugate();
            case 2:
                return z[startIdx[1] + getRecIdx(b,a)].conjugate();
            }
        case 2:
            switch(bt) {
            case 0:
                return z[startIdx[1] + getRecIdx(a,b)];
            case 1:
                return z[startIdx[1] + getRecIdx(a,b)].conjugate();
            case 2:
                return z[startIdx[4] + getTriIdx(a,b)];
            }
        }
        
        // This should never happen
        System.err.println("Invalid superscripts passed to cumulantDouble: " + al + " " + bt);
        
        return null;
    }
    
    private Complex cumulant(int al, int bt, int gm, int a, int b, int c, Complex[] z) {
        Complex ans = cumulantDouble(al, bt, a, b, z).multiply(cumulantSingle(gm, c, z));
        ans = ans.add(cumulantDouble(bt, gm, b, c, z).multiply(cumulantSingle(al, a, z)));
        ans = ans.add(cumulantDouble(al, gm, a, c, z).multiply(cumulantSingle(bt, b, z)));
        ans = ans.subtract(cumulantSingle(al, a, z).multiply(cumulantSingle(bt, b, z)).multiply(cumulantSingle(gm, c, z)).multiply(2.0));
        
        return ans;
    }

    @Override
    public void computeDerivatives(double t, Complex[] z, Complex[] zDot)
            throws MaxCountExceededException, DimensionMismatchException {
        int idx;
        
        // sigma_a^+
        for(int a = 0; a < n; ++a) {
            Complex t2 = new Complex(0);
            idx = startIdx[1] + a*(n-1);
            for(int b = 0; b < n; ++b) {
                if(a == b) continue;
                
                t2 = t2.add(coupling.getAlpha(a, b).conjugate().multiply(z[idx]));
                ++idx;
            }
            t2 = t2.multiply(0.5*gamma);
            
            zDot[a] = z[a].multiply(c1[a]).add(t2);
        }
        
        // sigma_a^z sigma_b^+
        idx = startIdx[1];
        for(int a = 0; a < n; ++a) {
            for(int b = 0; b < n; ++b) {
                if(a == b) continue;
                
                Complex t1 = new Complex(0);
                for(int j = 0; j < n; ++j) {
                    if(j == a || j == b) continue;
                    
                    t1 = t1.add(coupling.getAlpha(b, j).conjugate().multiply(cumulant(2, 2, 0, a, b, j, z)));
                }
                t1 = t1.multiply(0.5*gamma);
                
                Complex t2 = new Complex(0);
                for(int j = 0; j < n; ++j) {
                    if(j == a || j == b) continue;
                    
                    t2 = t2.add(coupling.getAlpha(a, j).multiply(cumulant(0, 0, 1, a, b, j, z)));
                }
                t2 = t2.multiply(-1.0*gamma);
                
                Complex t3 = new Complex(0);
                for(int j = 0; j < n; ++j) {
                    if(j == a || j == b) continue;
                    
                    t3 = t3.add(coupling.getAlpha(a, j).conjugate().multiply(cumulant(1, 0, 0, a, b, j, z)));
                }
                t3 = t3.multiply(-1.0*gamma);
                
                zDot[idx] = z[idx].multiply(c2[b]).add(z[b].multiply(-diff_gw)).add(coupling.getAlpha(a, b).multiply(z[a]).multiply(-0.5*gamma));
                zDot[idx] = zDot[idx].add(z[n+getRecIdx(b,a)].multiply(coupling.getF(a, b)).multiply(-gamma));
                zDot[idx] = zDot[idx].add(t1).add(t2).add(t3);
                ++idx;
            }
        }
        
        // sigma_a^z
        for(int a = 0; a < n; ++a) {
            Complex t1 = new Complex(0);
            for(int b = 0; b < n; ++b) {
                if(a == b) continue;
                
                t1 = t1.add(cumulantDouble(0, 1, a, b, z).subtract(cumulantDouble(1, 0, a, b, z)).multiply(coupling.getG(a, b)));
            }
            t1 = t1.multiply(migamma);
            
            Complex t2 = new Complex(0);
            for(int b = 0; b < n; ++b) {
                if(a == b) continue;
                
                t2 = t2.add(cumulantDouble(0, 1, a, b, z).add(cumulantDouble(1, 0, a, b, z)).multiply(coupling.getF(a, b)));
            }
            t2 = t2.multiply(-gamma);
            
            zDot[startIdx[2] + a] = z[startIdx[2] + a].multiply(-(gamma + w)).add(w - gamma).add(t1).add(t2);
        }
        
        // sigma_a^+ sigma_b^-
        idx = startIdx[3];
        for(int a = 0; a < n; ++a) {
            for(int b = a+1; b < n; ++b) {
                Complex t1 = new Complex(0);
                for(int j = 0; j < n; ++j) {
                    if(j == a || j == b) continue;
                    
                    t1 = t1.add(coupling.getAlpha(a, j).conjugate().multiply(cumulant(2, 1, 0, a, b, j, z)));
                }
                t1 = t1.multiply(0.5*gamma);
                
                Complex t2 = new Complex(0);
                for(int j = 0; j < n; ++j) {
                    if(j == a || j == b) continue;
                    
                    t2 = t2.add(coupling.getAlpha(b, j).multiply(cumulant(2, 0, 1, b, a, j, z)));
                }
                t2 = t2.multiply(0.5*gamma);
                
                zDot[idx] = z[idx].multiply(c3[a][b]).add(migamma.multiply(0.25).multiply(coupling.getG(a,b)).multiply(z[startIdx[2] + a].subtract(z[startIdx[2] + b])));
                zDot[idx] = zDot[idx].add(cumulantDouble(2, 2, a, b, z).add((z[startIdx[2] + a].add(z[startIdx[2] + b])).multiply(0.5)).multiply(0.5*gamma*coupling.getF(a,b)));
                zDot[idx] = zDot[idx].add(t1).add(t2);
                ++idx;
            }
        }
        
        // sigma_a^z sigma_b^z
        idx = startIdx[4];
        for(int a = 0; a < n; ++a) {
            for(int b = a+1; b < n; ++b) {
                Complex t1 = new Complex(0);
                for(int j = 0; j < n; ++j) {
                    if(j == a || j == b) continue;
                    
                    t1 = t1.add(cumulant(0, 2, 1, a, b, j, z).subtract(cumulant(1, 2, 0, a, b, j, z)).multiply(coupling.getG(a, j)));
                }
                t1 = t1.multiply(migamma);
                
                Complex t2 = new Complex(0);
                for(int j = 0; j < n; ++j) {
                    if(j == a || j == b) continue;
                    
                    t2 = t2.add(cumulant(0, 2, 1, b, a, j, z).subtract(cumulant(1, 2, 0, b, a, j, z)).multiply(coupling.getG(b, j)));
                }
                t2 = t2.multiply(migamma);
                
                Complex t3 = new Complex(0);
                for(int j = 0; j < n; ++j) {
                    if(j == a || j == b) continue;
                    
                    t3 = t3.add(cumulant(0, 2, 1, a, b, j, z).add(cumulant(1, 2, 0, a, b, j, z)).multiply(coupling.getF(a, j)));
                }
                t3 = t3.multiply(-gamma);
                
                Complex t4 = new Complex(0);
                for(int j = 0; j < n; ++j) {
                    if(j == a || j == b) continue;
                    
                    t4 = t4.add(cumulant(0, 2, 1, b, a, j, z).add(cumulant(1, 2, 0, b, a, j, z)).multiply(coupling.getF(b, j)));
                }
                t4 = t4.multiply(-gamma);
                
                zDot[idx] = z[startIdx[2] + a].add(z[startIdx[2] + b]).multiply(w - gamma);
                zDot[idx] = zDot[idx].subtract(cumulantDouble(2, 2, a, b, z).multiply(2*(gamma + w)));
                zDot[idx] = zDot[idx].add(cumulantDouble(0, 1, a, b, z).add(cumulantDouble(1, 0, a, b, z)).multiply(2*gamma*coupling.getF(a, b)));
                zDot[idx] = zDot[idx].add(t1).add(t2).add(t3).add(t4);
                ++idx;
            }
        }
        
        // sigma_a^+ sigma_b^+
        idx = startIdx[5];
        for(int a = 0; a < n; ++a) {
            for(int b = a+1; b < n; ++b) {
                Complex t1 = new Complex(0);
                for(int j = 0; j < n; ++j) {
                    if(j == a || j == b) continue;
                    
                    t1 = t1.add(cumulant(2, 0, 0, a, b, j, z).multiply(coupling.getAlpha(a, j).conjugate()));
                }
                t1 = t1.multiply(0.5*gamma);
                
                Complex t2 = new Complex(0);
                for(int j = 0; j < n; ++j) {
                    if(j == a || j == b) continue;
                    
                    t2 = t2.add(cumulant(2, 0, 0, b, a, j, z).multiply(coupling.getAlpha(b, j).conjugate()));
                }
                t2 = t2.multiply(0.5*gamma);
                
                zDot[idx] = z[idx].multiply(c4[a][b]).add(t1).add(t2);
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
