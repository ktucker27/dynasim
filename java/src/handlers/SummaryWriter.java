package handlers;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.Iterator;
import java.util.TreeMap;

import ode.CumulantParams;
import eval.SystemEval;

public class SummaryWriter {
    
    private class SummaryVals {
        double orderParam;
        double avgSigmaz;

        public SummaryVals() {
            orderParam = 0.0;
            avgSigmaz = 0.0;
        }
    }
    
    TreeMap<CumulantParams, SummaryVals> myVals;
    PrintWriter myWriter;
    boolean myLiveUpdate;
    boolean myHeaderOut;
    
    public SummaryWriter(String filepath) throws FileNotFoundException {
        myVals = new TreeMap<CumulantParams, SummaryWriter.SummaryVals>();
        myWriter = new PrintWriter(filepath);
        myLiveUpdate = false;
        myHeaderOut = false;
    }
    
    public void setLiveUpdate(boolean liveUpdate) {
        myLiveUpdate = liveUpdate;
    }
    
    public void addVals(CumulantParams params, SystemEval eval, double[] y) {
        SummaryVals vals = new SummaryVals();
        vals.orderParam = eval.getOrderParam(y);
        vals.avgSigmaz = eval.getAvgSigmaz(y);
        myVals.put(params, vals);
        
        if(myLiveUpdate) {
            writeVals(params, vals);
        }
    }
    
    public void addVals(CumulantParams params, DataRecorder recorder) {
        SummaryVals vals = new SummaryVals();
        vals.orderParam = recorder.getMeanOrderParam();
        vals.avgSigmaz = recorder.getMeanAvgZs();
        myVals.put(params, vals);
        
        if(myLiveUpdate) {
            writeVals(params, vals);
        }
    }
    
    public void close() {
        myWriter.close();
    }
    
    public void writeAllToFile() {
        // Write the header
        if(!myHeaderOut) {
            myWriter.write(CumulantParams.getHeader() + ", " + getHeader() + "\n");
            myHeaderOut = true;
        }
        
        // Write the values
        Iterator<CumulantParams> iter = myVals.keySet().iterator();
        while(iter.hasNext()) {
            CumulantParams params = iter.next();
            SummaryVals vals = myVals.get(params);
            myWriter.write(params.getLine() + ", ");
            myWriter.write(vals.orderParam + ", ");
            myWriter.write(vals.avgSigmaz + "\n");
        }
        
        myWriter.close();
    }
    
    private void writeVals(CumulantParams params, SummaryVals vals) {
        if(!myHeaderOut) {
            myWriter.write(CumulantParams.getHeader() + ", " + getHeader() + "\n");
            myHeaderOut = true;
        }
        
        myWriter.write(params.getLine() + ", ");
        myWriter.write(vals.orderParam + ", ");
        myWriter.write(vals.avgSigmaz + "\n");
        
        myWriter.flush();
    }
    
    private String getHeader() {
        return "order_param, avg_sigmaz";
    }
}
