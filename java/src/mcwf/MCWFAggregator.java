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
    private ExpectedSpinValues[] myEvs;
    
    private final double TIME_TOL = 1e-10;
    
    public MCWFAggregator(int numTimes, double timeDelta) {
        myNumTrajectories = 0;
        myNumTimes = numTimes;
        myTimeDelta = timeDelta;
        
        myEvs = new ExpectedSpinValues[numTimes];
        for(int i = 0; i < numTimes; ++i) {
            myEvs[i] = new ExpectedSpinValues();
        }
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
    
    public ExpectedSpinValues getEvs(int idx) {
        return myEvs[idx];
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
                myEvs[timeIdx].addEq(traj.getEvs(timeIdx));
            }
        }
    }
}
