package mcwf;

import java.util.ArrayList;

public class TrajectoryStats {
    ArrayList<Double> myStopTimes;
    ArrayList<Jump> myJumps;

    private static final double TIME_TOL = 1.0e-10;
    
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
    
    public TrajectoryStats() {
        myStopTimes = new ArrayList<Double>();
        myJumps = new ArrayList<Jump>();
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
    
    /**
     * Performs a number of consistency checks on the stats for this trajectory
     * 
     * @param msgs list of error messages
     * @return True if stats are consistent
     */
    public boolean validate(ArrayList<String> msgs) {
        boolean valid = true;
        if(myStopTimes.size() != myJumps.size()) {
            msgs.add(String.format("Number of stop times %d not equal to number of jumps %d", myStopTimes.size(), myJumps.size()));
            valid = false;
        } else {
            for(int i = 0; i < myStopTimes.size(); ++i) {
                if(Math.abs(myStopTimes.get(i) - getJumpTime(i)) > TIME_TOL) {
                    msgs.add(String.format("Stop time %d: %g different from jump time: %g", i, myStopTimes.get(i), getJumpTime(i)));
                    valid = false;
                }
            }
        }
        
        return valid;
    }
}
