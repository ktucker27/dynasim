package mcwf;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import utils.ExpectedSpinValues;
import utils.HusimiDist;

/**
 * Class that computes and stores aggregate expected values by time coming out of
 * the MCWFIntegrator
 */
public class MCWFAggregator {
    // Simulation parameters
    private int myNumTrajectories;
    private int myNumTimes;
    private double myTimeDelta;
    
    // Output data indexed by time
    private ExpectedSpinValues[] mySumEvs;
    private ExpectedSpinValues[] mySumSqEvs;
    
    // Optional aggregated stats
    private int[][] myNumJumps;
    private int[] mySumJs;
    HashMap<Double, HusimiDist> myHusimis;
    
    private final double TIME_TOL = 1e-10;
    
    public MCWFAggregator(int numTimes, double timeDelta) {
        myNumTrajectories = 0;
        myNumTimes = numTimes;
        myTimeDelta = timeDelta;
        
        mySumEvs = new ExpectedSpinValues[numTimes];
        mySumSqEvs = new ExpectedSpinValues[numTimes];
        for(int i = 0; i < numTimes; ++i) {
            mySumEvs[i] = new ExpectedSpinValues();
            mySumSqEvs[i] = new ExpectedSpinValues();
        }
        
        myNumJumps = null;
        mySumJs = null;
        myHusimis = null;
    }
    
    public int getNumTrajectories() {
        return myNumTrajectories;
    }
    
    public int getNumTimes() {
        return myNumTimes;
    }
    
    public double getTime(int idx) {
        return idx*myTimeDelta;
    }
    
    public ExpectedSpinValues getSumEvs(int idx) {
        return mySumEvs[idx];
    }
    
    public ExpectedSpinValues getSumSqEvs(int idx) {
        return mySumSqEvs[idx];
    }
    
    public int[] getJumps(int idx) {
        return myNumJumps[idx];
    }
    
    public int getSumJs(int idx) {
        return mySumJs[idx];
    }
    
    public Map<Double, HusimiDist> getHusimis() {
        return myHusimis;
    }
    
    public void aggregate(QuantumTrajectory[] trajectories, int startIdx, int endIdx) {
        for(int idx = startIdx; idx <= endIdx; ++idx) {
            ++myNumTrajectories;
            
            QuantumTrajectory traj = trajectories[idx];
            if(myNumTimes != traj.getNumEvTimes() ||
               Math.abs(myTimeDelta - traj.getEvTimeDelta()) > TIME_TOL) {
                throw new UnsupportedOperationException("MCWFAggregator: Trajectory does not align with time grid " + myNumTimes + " " + traj.getNumEvTimes() + " " + myTimeDelta + " " + traj.getEvTimeDelta());
            }
            
            for(int timeIdx = 0; timeIdx < myNumTimes; ++timeIdx) {
                if(Math.abs(traj.getEvs(timeIdx).getTime() - timeIdx*myTimeDelta) > TIME_TOL) {
                    throw new UnsupportedOperationException("MCWFAggregator: Expected value time " + traj.getEvs(timeIdx).getTime() + " does not match expected time " + timeIdx*myTimeDelta);
                }
                mySumEvs[timeIdx].addEq(traj.getEvs(timeIdx));
                mySumSqEvs[timeIdx].addSqEq(traj.getEvs(timeIdx));
            }
            
            // Aggregate optional stats
            if(traj.getStats() != null) {
                aggregateStats(traj);
            }
        }
    }
    
    private void aggregateStats(QuantumTrajectory traj) {
        // Initialize the appropriate structures if this is the first trajectory seen
        TrajectoryStats stats = traj.getStats();
        boolean first = false;
        if(myNumJumps == null) {
            first = true;
            
            myNumJumps = new int[myNumTimes][traj.getNumOutcomes()-1];
            mySumJs = new int[myNumTimes];
            
            for(int timeIdx = 0; timeIdx < myNumTimes; ++timeIdx) {
                for(int jumpIdx = 0; jumpIdx < traj.getNumOutcomes() - 1; ++jumpIdx) {
                    myNumJumps[timeIdx][jumpIdx] = 0;
                }
                mySumJs[timeIdx] = 0;
            }
            
            myHusimis = new HashMap<Double, HusimiDist>();
            Map<Double, HusimiDist> husimis = stats.getHusimis();
            Iterator<Map.Entry<Double, HusimiDist>> iter = husimis.entrySet().iterator();
            while(iter.hasNext()) {
                Map.Entry<Double, HusimiDist> entry = iter.next();
                myHusimis.put(entry.getKey(), new HusimiDist(entry.getValue()));
            }
        }
        
        // Aggregate number of jumps
        int numJumps = stats.getNumJumps();
        for(int jumpIdx = 0; jumpIdx < numJumps; ++jumpIdx) {
            double jumpTime = stats.getJumpTime(jumpIdx);
            int outcome = stats.getJumpOutcome(jumpIdx);
            
            int jumpTimeIdx = (int)(jumpTime/myTimeDelta);
            ++myNumJumps[jumpTimeIdx][outcome];
        }
        
        // Aggregate J values
        for(int timeIdx = 0; timeIdx < myNumTimes; ++timeIdx) {
            mySumJs[timeIdx] += stats.getJ(timeIdx);
        }
        
        // Aggregate Husimi distributions
        if(!first) {
            Map<Double, HusimiDist> husimis = stats.getHusimis();
            if(myHusimis.size() != husimis.size()) {
                throw new UnsupportedOperationException("Found trajectory with inconsistent number of Husimi distributions");
            }
            
            Iterator<Map.Entry<Double, HusimiDist>> iter = husimis.entrySet().iterator();
            while(iter.hasNext()) {
                Map.Entry<Double, HusimiDist> entry = iter.next();
                HusimiDist dist = myHusimis.get(entry.getKey());
                if(dist == null) {
                    throw new UnsupportedOperationException("Found Husimi distribution at inconsistent time");
                }
                dist.addEq(entry.getValue());
            }
        }
    }
}
