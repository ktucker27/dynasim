package ode;

public class RPAEval {
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

    public int getDimension() {
        return 3*n + 3*n*(n-1)/2 + 3*n*(n-1);
    }
    
    public int getN() {
        return n;
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
        if(al == 2 && bt != 2 && gm != 2) {
            return rpaSum(al, bt, gm, a, b, 0);
        } else if(al != 2 && bt == 2 && gm != 2) {
            return rpaSum(al, bt, gm, a, b, 1);
        } else if(al != 2 && bt != 2 && gm == 2) {
            return rpaSum(al, bt, gm, a, b, 2);
        } else {
            return cumulantSum(al, bt, gm, a, b);
        }
    }
    
    public double compCorr() {
        double sum = 0.0;
        for(int a = 0; a < n; ++a) {
            sum += rowsums[2][2][a] + rowsums[1][1][a];
        }
        
        return 0.25*sum/(double)(n*(n-1));
    }
}
