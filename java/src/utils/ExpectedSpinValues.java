package utils;

public class ExpectedSpinValues {
    private double myTime;
    private double[] myEs;
    private DynaComplex[][] myEss;
    
    /**
     * Initializes all expected values to zero
     */
    public ExpectedSpinValues() {
        myTime = 0.0;
        
        myEs = new double[3];
        
        myEss = new DynaComplex[3][3];
        for(int i = 0; i < 3; ++i) {
            myEs[i] = 0.0;
            for(int j = 0; j < 3; ++j) {
                myEss[i][j] = new DynaComplex(0,0);
            }
        }
    }
    
    public void setTime(double time) {
        myTime = time;
    }
    
    public double getTime() {
        return myTime;
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
