package utils;

import handlers.SummaryWriter;
import integrator.IntegratorReducer;

public class SynchReducer implements IntegratorReducer {

    private SummaryWriter myWriter;
    
    public SynchReducer(SummaryWriter writer) {
        myWriter = writer;
    }
    
    @Override
    public void reduce(ODESolution soln) {
        if(!(soln instanceof SynchSolution)) {
            throw new UnsupportedOperationException("SynchReducer can only reduce solutions of type SynchSolution");
        }
        
        SynchSolution synchSoln = (SynchSolution)soln;
        myWriter.addVals(synchSoln.getParams(), synchSoln.getRecorder());
    }

}
