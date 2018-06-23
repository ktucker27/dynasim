package utils;

import eval.SystemEval;
import handlers.DataRecorder;
import ode.SystemParams;

public class SynchSolution implements ODESolution {

    private SystemParams myParams;
    private DataRecorder myRecorder;
    private SystemEval myEval;
    private double[] y;

    public SynchSolution(SystemParams params, DataRecorder recorder, SystemEval eval) {
        super();
        myParams = params;
        myRecorder = recorder;
        myEval = eval;
        y = null;
    }
    
    public SystemParams getParams() {
        return myParams;
    }
    
    public DataRecorder getRecorder() {
        return myRecorder;
    }
    
    public SystemEval getEval() {
        return myEval;
    }

    public double[] getSolution() {
        return y;
    }

    @Override
    public void setSolution(double[] y) {
        this.y = y;
    }

}
