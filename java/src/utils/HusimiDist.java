package utils;

import org.apache.commons.math3.util.CombinatoricsUtils;

public class HusimiDist {

    private int myNumTheta;
    private int myNumPhi;
    private double myDeltaTheta;
    private double myDeltaPhi;
    private boolean myIsZero;
    
    private double[][] myVals;
    
    public HusimiDist(int numTheta, int numPhi) {
        double deltaTheta = Math.PI/((double)(numTheta-1));
        double deltaPhi = 2.0*Math.PI/((double)(numPhi-1));
        
        myNumTheta = numTheta;
        myNumPhi = numPhi;
        myDeltaTheta = deltaTheta;
        myDeltaPhi = deltaPhi;
        myIsZero = true;
        
        myVals = new double[numTheta][numPhi];
        for(int i = 0; i < numTheta; ++i) {
            for(int j = 0; j < numPhi; ++j) {
                myVals[i][j] = 0.0;
            }
        }
    }
    
    public HusimiDist(HusimiDist rhs) {
        this.myNumTheta = rhs.myNumTheta;
        this.myNumPhi = rhs.myNumPhi;
        this.myDeltaTheta = rhs.myDeltaTheta;
        this.myDeltaPhi = rhs.myDeltaPhi;
        this.myIsZero = rhs.myIsZero;
        
        myVals = new double[myNumTheta][myNumPhi];
        for(int i = 0; i < myNumTheta; ++i) {
            for(int j = 0; j < myNumPhi; ++j) {
                myVals[i][j] = rhs.myVals[i][j];
            }
        }
    }
    
    public boolean isZero() {
        return myIsZero;
    }
    
    public int getNumTheta() {
        return myNumTheta;
    }
    
    public int getNumPhi() {
        return myNumPhi;
    }
    
    public final double[][] getVals() {
        return myVals;
    }
    
    public void calc(DynaComplex[] state, int jval) {
        int n = state.length - 1;
        
        // TODO - Is the distribution zero for less than maximal total spin?
        if(jval != n/2) {
            myIsZero = true;
            return;
        }
        
        DynaComplex ip = new DynaComplex(0,0);
        DynaComplex t1 = new DynaComplex(0,0);
        DynaComplex t2 = new DynaComplex(0,0);
        
        for(int thetaIdx = 0; thetaIdx < myNumTheta; ++thetaIdx) {
            double theta = thetaIdx*myDeltaTheta;
            double ct2 = Math.cos(0.5*theta);
            double st2 = Math.sin(0.5*theta);
            for(int phiIdx = 0; phiIdx < myNumPhi; ++phiIdx) {
                double phi = phiIdx*myDeltaPhi;
                
                // Compute the inner product between the state and the CSS
                ip.set(0,0);
                for(int i = 0; i < state.length; ++i) {
                    double ct2pow = Math.pow(ct2, n - i);
                    double st2pow = Math.pow(st2, i);
                    
                    t1.set(Math.cos(i*phi), -Math.sin(i*phi));
                    CombinatoricsUtils.checkBinomial(n, i);
                    t2.set(ct2pow, 0).multiply(st2pow).multiply(t1).multiply(Math.sqrt(CombinatoricsUtils.binomialCoefficientDouble(n, i)));
                    ip.add(t1.set(state[i]).multiply(t2));
                }
                myVals[thetaIdx][phiIdx] = ip.modSq();
            }
        }
        
        myIsZero = false;
    }
    
    public void addEq(HusimiDist rhs) {
        if(myNumTheta != rhs.myNumTheta || myNumPhi != rhs.myNumPhi) {
            throw new UnsupportedOperationException("Trying to add incompatibly sized Husimi distributions");
        }
        
        if(rhs.myIsZero) return;
        
        for(int i = 0; i < myNumTheta; ++i) {
            for(int j = 0; j < myNumPhi; ++j) {
                myVals[i][j] += rhs.myVals[i][j];
            }
        }
        myIsZero = false;
    }
}
