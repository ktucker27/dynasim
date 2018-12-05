package mcwf;

import utils.ExpectedSpinValues;

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
    
    public void aggregate(QuantumTrajectory[] trajectories, int startIdx, int endIdx) {
        for(int idx = startIdx; idx <= endIdx; ++idx) {
            ++myNumTrajectories;
            
            QuantumTrajectory traj = trajectories[idx];
            if(myNumTimes != traj.getNumEvTimes() ||
               myTimeDelta != traj.getEvTimeDelta()) {
                throw new UnsupportedOperationException("MCWFAggregator: Trajectory does not align with time grid");
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
        if(myNumJumps == null) {
            myNumJumps = new int[myNumTimes][traj.getNumOutcomes()-1];
            mySumJs = new int[myNumTimes];
            
            for(int timeIdx = 0; timeIdx < myNumTimes; ++timeIdx) {
                for(int jumpIdx = 0; jumpIdx < traj.getNumOutcomes() - 1; ++jumpIdx) {
                    myNumJumps[timeIdx][jumpIdx] = 0;
                }
                mySumJs[timeIdx] = 0;
            }
        }
        
        TrajectoryStats stats = traj.getStats();
        int numJumps = stats.getNumJumps();
        for(int jumpIdx = 0; jumpIdx < numJumps; ++jumpIdx) {
            double jumpTime = stats.getJumpTime(jumpIdx);
            int outcome = stats.getJumpOutcome(jumpIdx);
            
            int jumpTimeIdx = (int)(jumpTime/myTimeDelta);
            ++myNumJumps[jumpTimeIdx][outcome];
        }
        
        for(int timeIdx = 0; timeIdx < myNumTimes; ++timeIdx) {
            mySumJs[timeIdx] += stats.getJ(timeIdx);
        }
    }
}
