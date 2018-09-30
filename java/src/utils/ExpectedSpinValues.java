package utils;

public class ExpectedSpinValues {
    private double[] myEs;
    private DynaComplex[][] myEss;
    
    /**
     * Initializes all expected values to zero
     */
    public ExpectedSpinValues() {
        myEs = new double[3];
        
        myEss = new DynaComplex[3][3];
        for(int i = 0; i < 3; ++i) {
            myEs[i] = 0.0;
            for(int j = 0; j < 3; ++j) {
                myEss[i][j] = new DynaComplex(0,0);
            }
        }
    }

    public void setEs(int i, double val) {
        myEs[i] = val;
    }

    public double getEs(int i) {
        return myEs[i];
    }
    
    public DynaComplex getEss(int i, int j) {
        return myEss[i][j];
    }

    public void addEq(ExpectedSpinValues rhs) {
        for(int i = 0; i < 3; ++i) {
            myEs[i] += rhs.myEs[i];
            for(int j = 0; j < 3; ++j) {
                myEss[i][j].add(rhs.myEss[i][j]);
            }
        }
    }
}
