package mcwf;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import utils.DynaComplex;
import utils.HusimiDist;

public class TrajectoryStats {
    ArrayList<Double> myStopTimes;
    ArrayList<Jump> myJumps;
    HashMap<Double, HusimiDist> myHusimis;
    int[] myJs;

    private static final double TIME_TOL = 1.0e-10;
    private static final int NUM_THETA = 180;
    private static final int NUM_PHI = 360;
    
    private class Jump {
        private double myTime;
        private int myOutcome;
        
        public Jump(double time, int outcome) {
            myTime = time;
            myOutcome = outcome;
        }
        
        public double getTime() {
            return myTime;
        }
        
        public int getOutcome() {
            return myOutcome;
        }
    }
    
    public TrajectoryStats(int numTimes) {
        myStopTimes = new ArrayList<Double>();
        myJumps = new ArrayList<Jump>();
        
        myHusimis = new HashMap<Double, HusimiDist>();
        myHusimis.put(1.0, new HusimiDist(NUM_THETA, NUM_PHI));
        
        myJs = new int[numTimes];
        for(int i = 0; i < numTimes; ++i) {
            myJs[i] = 0;
        }
    }
    
    public void recordStopTime(double time) {
        myStopTimes.add(time);
    }
    
    public int getNumStopTimes() {
        return myStopTimes.size();
    }
    
    public double getStopTime(int idx) {
        return myStopTimes.get(idx);
    }
    
    public void recordJump(double time, int outcome) {
        myJumps.add(new Jump(time, outcome));
    }
    
    public int getNumJumps() {
        return myJumps.size();
    }
    
    public double getJumpTime(int idx) {
        return myJumps.get(idx).getTime();
    }
    
    public int getJumpOutcome(int idx) {
        return myJumps.get(idx).getOutcome();
    }
    
    public int getJ(int idx) {
        return myJs[idx];
    }
    
    public Map<Double, HusimiDist> getHusimis() {
        return myHusimis;
    }
    
    public void recordState(int timeIdx, double time, DynaComplex[] state, int jval) {
        myJs[timeIdx] = jval;
        
        if(myHusimis.size() > 0) {
            Iterator<Map.Entry<Double, HusimiDist>> iter = myHusimis.entrySet().iterator();
            while(iter.hasNext()) {
                Map.Entry<Double, HusimiDist> entry = iter.next();
                if(Math.abs(entry.getKey() - time) < TIME_TOL) {
                    entry.getValue().calc(state, jval);
                }
            }
        }
    }
    
    /**
     * Performs a number of consistency checks on the stats for this trajectory
     * 
     * @param msgs list of error messages
     * @return True if stats are consistent
     */
    public boolean validate(ArrayList<String> msgs) {
        boolean valid = true;
        if(myStopTimes.size() - 1 != myJumps.size()) {
            msgs.add(String.format("Number of stop times %d not equal to number of jumps %d", myStopTimes.size(), myJumps.size()));
            valid = false;
        } else {
            for(int i = 0; i < myStopTimes.size() - 1; ++i) {
                if(Math.abs(myStopTimes.get(i) - getJumpTime(i)) > TIME_TOL) {
                    msgs.add(String.format("Stop time %d: %g different from jump time: %g", i, myStopTimes.get(i), getJumpTime(i)));
                    valid = false;
                }
            }
        }
        
        return valid;
    }
}
