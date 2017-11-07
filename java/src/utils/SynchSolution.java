package utils;

import eval.SystemEval;
import handlers.DataRecorder;
import ode.CumulantParams;

public class SynchSolution implements ODESolution {

    private CumulantParams myParams;
    private DataRecorder myRecorder;
    private SystemEval myEval;
    private double[] y;

    public SynchSolution(CumulantParams params, DataRecorder recorder, SystemEval eval) {
        super();
        myParams = params;
        myRecorder = recorder;
        myEval = eval;
        y = null;
    }
    
    public CumulantParams getParams() {
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
