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
                throw new UnsupportedOperationException("Trajectory does not align with time grid");
            }
            
            for(int timeIdx = 0; timeIdx < myNumTimes; ++timeIdx) {
                myEvs[timeIdx].addEq(traj.getEvs(timeIdx));
            }
        }
    }
}
