package utils;

import java.io.FileNotFoundException;

import handlers.SummaryWriter;
import integrator.IntegratorReducer;

public class SynchReducer implements IntegratorReducer {

    private SummaryWriter myWriter;
    private String myFinalDir;
    
    public SynchReducer() {
        myWriter = null;
        myFinalDir = new String();
    }
    
    public void writeSummaries(SummaryWriter writer) {
        myWriter = writer;
    }
    
    public void writeFinal(String finalDir) {
        myFinalDir = finalDir;
    }
    
    @Override
    public void reduce(ODESolution soln) {
        if(!(soln instanceof SynchSolution)) {
            throw new UnsupportedOperationException("SynchReducer can only reduce solutions of type SynchSolution");
        }
        
        SynchSolution synchSoln = (SynchSolution)soln;
        
        if(myWriter != null) {
            myWriter.addVals(synchSoln.getParams(), synchSoln.getRecorder());
        }
        
        if(!myFinalDir.isEmpty()) {
            String filepath = myFinalDir + "/final_" + synchSoln.getParams().getFilename();
            try {
                SynchUtils.writeToFile(filepath, synchSoln.getSolution());
            } catch(FileNotFoundException ex) {
                System.err.println("ERROR - Could not write solution file " + filepath);
            }
        }
    }
}
