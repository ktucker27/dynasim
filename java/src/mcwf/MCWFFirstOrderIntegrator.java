package mcwf;

public class MCWFFirstOrderIntegrator implements Runnable {
    // Simulation parameters
    private QuantumTrajectory myTrajectory;
    
    public MCWFFirstOrderIntegrator(QuantumTrajectory trajectory) {
        myTrajectory = trajectory;
    }
    
    @Override
    public void run() {
        // Populate expected values at the start time
        myTrajectory.calcEvs(0);
        
        double[] cdf = new double[myTrajectory.getNumOutcomes()-1];
        
        for(int timeIdx = 1; timeIdx < myTrajectory.getNumTimes(); ++timeIdx) {
            // Determine the jump
            myTrajectory.getCdf(cdf);
            int outcome = myTrajectory.rollDice(cdf);
            
            myTrajectory.performJump(outcome);
            
            // Calculate the expected values
            if(timeIdx % myTrajectory.getNumStepsPerEv() == 0) {
                myTrajectory.calcEvs(timeIdx*myTrajectory.getTimeDelta());
            }
        }
    }
}
