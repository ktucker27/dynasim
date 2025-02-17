package eval;

import java.util.Random;

import utils.DynaComplex;


public class RPAEval implements SystemEval {
    int n;
    
    double[] myVals;
    
    /**
     * The following convention is used for all Pauli operator types in this class:
     * 0 = x, 1 = y, 2 = z  
     */
    
    /**
     * vsums[i] = the sum of all <s_a^i>
     */
    double[] vsums;
    
    /**
     * rowsums[i][j][a] = sum_b <s_a^i s_b^j>
     */
    double[][][] rowsums;
    
    /**
     * Stores the start indices for each variable establishing the following order:
     * 0: <s_a^x>
     * 1: <s_a^y>
     * 2: <s_a^z>
     * 3: <s_a^x s_b^x>
     * 4: <s_a^y s_b^y>
     * 5: <s_a^z s_b^z>
     * 6: <s_a^x s_b^y>
     * 7: <s_a^x s_b^z>
     * 8: <s_a^y s_b^z>
     */
    int[] startIdx;
    
    public RPAEval(int n) {
        this.n = n;
        
        init();
    }
    
    private void init() {
        startIdx = new int[9];
        startIdx[0] = 0;
        startIdx[1] = n;
        startIdx[2] = startIdx[1] + n;
        startIdx[3] = startIdx[2] + n;
        startIdx[4] = startIdx[3] + n*(n-1)/2;
        startIdx[5] = startIdx[4] + n*(n-1)/2;
        startIdx[6] = startIdx[5] + n*(n-1)/2;
        startIdx[7] = startIdx[6] + n*(n-1);
        startIdx[8] = startIdx[7] + n*(n-1);
        
        vsums = new double[3];
        rowsums = new double[3][3][n];
    }
    
    public void setVals(double[] y) {
        this.myVals = y;
        
        // Start by computing the sums that we need.  This is O(n^2)

        // Zero out row and column sums
        for(int i = 0; i < 3; ++i) {
            vsums[i] = 0;
            for(int j = 0; j < 3; ++j) {
                for(int k = 0; k < n; ++k) {
                    rowsums[i][j][k] = 0;
                }
            }
        }

        for(int a = 0; a < n; ++a) {
            for(int i = 0; i < 3; ++i) {
                vsums[i] += y[startIdx[i] + a];
            }

            for(int b = 0; b < n; ++b) {
                if(a == b) continue;
                
                for(int i = 0; i < 3; ++i) {
                    for(int j = 0; j < 3; ++j) {
                        rowsums[i][j][a] += getDouble(i, j, a, b, y);
                    }
                }
            }
        }
    }
    
    public double[][][] getRowSums() {
        return rowsums;
    }
    
    public int[] getStartIdx() {
        return startIdx;
    }

    @Override
    public int getDimension() {
        return 3*n + 3*n*(n-1)/2 + 3*n*(n-1);
    }
    
    @Override
    public int getRealDimension() {
        return getDimension();
    }
    
    @Override
    public int getN() {
        return n;
    }
    
    @Override
    public double getOrderParam(double[] y) {
        setVals(y);
        return compCorr();
    }
    
    public double getAvgSigmax(double[] y) {
        double sum = 0.0;
        for(int a = 0; a < n; ++a) {
            sum += getSingle(2, a, y);
        }
        
        return sum/(double)n;
    }
    
    public double getAvgSigmay(double[] y) {
        double sum = 0.0;
        for(int a = 0; a < n; ++a) {
            sum += getSingle(1, a, y);
        }
        
        return sum/(double)n;
    }
    
    @Override
    public double getAvgSigmaz(double[] y) {
        double sum = 0.0;
        for(int a = 0; a < n; ++a) {
            sum -= getSingle(0, a, y);
        }
        
        return sum/(double)n;
    }
    
    @Override
    public double getAvgSigmazz(double[] y) {
        double sum = 0.0;
        for(int a = 0; a < n; ++a) {
            for(int b = a+1; b < n; ++b) {
                sum += getSingle(0, a, y)*getSingle(0, b, y);
            }
        }
        
        return sum/(double)(n*(n-1)/2);
    }
    
    public double getAvg(int al, double[] y) {
        double sum = 0.0;
        for(int a = 0; a < n; ++a) {
            sum += getSingle(fromGlobal(al), a, y);
        }
        
        return (al == 2 ? -1.0 : 1.0)*sum/(double)n;
    }
    
    public double getAvg(int al, int bt, double[] y) {
        setVals(y);
        return getAvg(al, bt);
    }
    
    public double getAvg(int al, int bt) {
        double amult = (al == 2 ? -1.0 : 1.0);
        double bmult = (bt == 2 ? -1.0 : 1.0);
        int lal = fromGlobal(al);
        int lbt = fromGlobal(bt);
        
        double sum = 0.0;
        for(int a = 0; a < n; ++a) {
            sum += rowsums[lal][lbt][a];
        }
        
        return amult*bmult*sum/(double)(n*(n-1));
    }
    
    private int fromGlobal(int al) {
        switch(al) {
        case 0:
            return 2;
        case 2:
            return 0;
        }
        
        return al;
    }
    
    @Override
    public void getBlochVectors(double[] y, double[] xs, double[] ys, double[] zs) {
        for(int a = 0; a < n; ++a) {
            xs[a] = getSingle(2, a, y);
            ys[a] = getSingle(1, a, y);
            zs[a] = -1.0*getSingle(0, a, y);
        }
    }
    
    @Override
    public void getFirstOrderCollectiveEvs(double[] y, double[] es) {
        throw new UnsupportedOperationException("RPAEval has not implemented getFirstOrderEvs");
    }
    
    @Override
    public void getSecondOrderCollectiveEvs(double[] y, DynaComplex[][] es) {
        throw new UnsupportedOperationException("RPAEval has not implemented getSecondOrderEvs");
    }
    
    @Override
    public void initSpinUpX(double[] y0) {
        initialize(y0, 0, 0, InitAngleType.CONST, InitAngleType.CONST);
    }
    
    @Override
    public void initialize(double[] y0, double zenith, double phase, InitAngleType zenithType, InitAngleType phaseType) {
        double zi = 0.0;
        double pi = 0.0;
        double szi = 0.0;
        double czi = 0.0;
        
        double eps = 1.0e-3;
        
        Random rg = new Random(1);
        for(int i = 0; i < n; ++i) {
            switch(zenithType) {
            case EQUAL_SPACING:
                zi = eps + i*(Math.PI - 2*eps)/(double)(n-1);
                break;
            case RANDOM:
                zi = rg.nextDouble()*Math.PI;
                break;
            case CONST:
                zi = zenith;
                break;
            }
            
            szi = Math.sin(zi);
            czi = Math.cos(zi);
            
            switch(phaseType) {
            case EQUAL_SPACING:
                pi = i*2*Math.PI/(double)n - Math.PI;
                break;
            case RANDOM:
                pi = rg.nextDouble()*2.0*Math.PI - Math.PI;
                break;
            case CONST:
                pi = phase;
                break;
            }
            
            y0[i] = -czi;
            y0[startIdx[1]+i] = szi*Math.sin(pi);
            y0[startIdx[2]+i] = szi*Math.cos(pi);
        }

        int triIdx = 0;
        int recIdx = 0;
        for(int a = 0; a < n; ++a) {
            for(int b = 0; b < n; ++b) {
                if(a == b) continue;
                
                if(b > a) {
                    y0[startIdx[3] + triIdx] = y0[a]*y0[b];
                    y0[startIdx[4] + triIdx] = y0[startIdx[1]+a]*y0[startIdx[1]+b];
                    y0[startIdx[5] + triIdx] = y0[startIdx[2]+a]*y0[startIdx[2]+b];
                    ++triIdx;
                }
                
                y0[startIdx[6] + recIdx] = y0[a]*y0[startIdx[1]+b];
                y0[startIdx[7] + recIdx] = y0[a]*y0[startIdx[2]+b];
                y0[startIdx[8] + recIdx] = y0[startIdx[1]+a]*y0[startIdx[2]+b];
                
                ++recIdx;
            }
        }
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
    
    /**
     * PRECONDITION: i < j
     * 
     * @param i First superscript
     * @param j Second superscript
     * @return Index offset for rectangular data array
     */
    private int getSuperTriIdx(int i, int j) {
        return 3*i - i*(i+1)/2 + j - i - 1;
    }
    
    /**
     * Does not require a call to setVals
     * 
     * @param al superscript
     * @param a particle index
     * @param y value array
     * @return <sigma_a^al> taken from y
     */
    public double getSingle(int al, int a, double[] y) {
        return y[startIdx[al] + a];
    }
    
    /**
     * Does not require a call to setVals
     * 
     * @param al first superscript
     * @param bt second superscript
     * @param a first particle index
     * @param b second particle index
     * @param y value array
     * @return <sigma_a^al sigma_b^bt> taken from y
     */
    public double getDouble(int al, int bt, int a, int b, double[] y) {
        if(al == bt) {
            return y[startIdx[3 + al] + getTriIdx(a,b)];
        }
        
        if(al < bt) {
            return y[startIdx[6 + getSuperTriIdx(al, bt)] + getRecIdx(a,b)];
        }
        
        return y[startIdx[6 + getSuperTriIdx(bt, al)] + getRecIdx(b,a)];
    }
    
    private double rpaSum(int al, int bt, int gm, int a, int b, int sidx) {
        switch(sidx) {
        case 0:
            return myVals[startIdx[al] + a]*(rowsums[bt][gm][b] - getDouble(bt, gm, b, a, myVals));
        case 1:
            return myVals[startIdx[bt] + b]*(rowsums[al][gm][a] - getDouble(al, gm, a, b, myVals));
        case 2:
            return (vsums[gm] - myVals[startIdx[gm] + a] - myVals[startIdx[gm] + b])*getDouble(al, bt, a, b, myVals);
        }
        
        // This should never happen
        System.err.println("Invalid sidx passed to rpaSum: " + sidx);
        
        return 0.0;
    }
    
    private double cumulantSum(int al, int bt, int gm, int a, int b) {
        double singleSum = vsums[gm] - myVals[startIdx[gm] + a] - myVals[startIdx[gm] + b];
        return (getDouble(al, bt, a, b, myVals)*singleSum
                + (rowsums[al][gm][a] - getDouble(al, gm, a, b, myVals))*myVals[startIdx[bt] + b]
                + (rowsums[bt][gm][b] - getDouble(bt, gm, b, a, myVals))*myVals[startIdx[al] + a]
                - 2*(myVals[startIdx[al] + a]*myVals[startIdx[bt] + b]*singleSum));
    }
    
    public double evalTriple(int al, int bt, int gm, int a, int b) {
        if(al == 0 && bt != 0 && gm != 0) {
            return rpaSum(al, bt, gm, a, b, 0);
        } else if(al != 0 && bt == 0 && gm != 0) {
            return rpaSum(al, bt, gm, a, b, 1);
        } else if(al != 0 && bt != 0 && gm == 0) {
            return rpaSum(al, bt, gm, a, b, 2);
        } else if(al == 0 && bt == 0 && gm == 0){
            return cumulantSum(al, bt, gm, a, b);
        }
        
        return 0.0;
    }
    
    public double compCorr() {
        double sum = 0.0;
        for(int a = 0; a < n; ++a) {
            sum += rowsums[2][2][a] + rowsums[1][1][a];
        }
        
        return 0.25*sum/(double)(n*(n-1));
    }
}
