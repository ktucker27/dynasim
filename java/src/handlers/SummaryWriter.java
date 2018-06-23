package handlers;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Iterator;
import java.util.TreeMap;

import eval.SystemEval;
import ode.SystemParams;

public class SummaryWriter {
    
    private class SummaryVals {
        double orderParam;
        double avgSigmaz;
        double avgSigmazz;

        public SummaryVals() {
            orderParam = 0.0;
            avgSigmaz = 0.0;
            avgSigmazz = 0.0;
        }
    }
    
    TreeMap<SystemParams, SummaryVals> myVals;
    PrintWriter myWriter;
    boolean myLiveUpdate;
    boolean myHeaderOut;
    
    public SummaryWriter(String filepath) throws IOException {
        myVals = new TreeMap<SystemParams, SummaryWriter.SummaryVals>();
        //myWriter = new PrintWriter(filepath);
        File file = new File(filepath);
        myHeaderOut = file.exists();
        FileWriter fw = new FileWriter(file, true);
        BufferedWriter bw = new BufferedWriter(fw);
        myWriter = new PrintWriter(bw);
        myLiveUpdate = false;
    }
    
    public void setLiveUpdate(boolean liveUpdate) {
        myLiveUpdate = liveUpdate;
    }
    
    public void addVals(SystemParams params, SystemEval eval, double[] y) {
        SummaryVals vals = new SummaryVals();
        vals.orderParam = eval.getOrderParam(y);
        vals.avgSigmaz = eval.getAvgSigmaz(y);
        myVals.put(params, vals);
        
        if(myLiveUpdate) {
            writeVals(params, vals);
        }
    }
    
    public void addVals(SystemParams params, DataRecorder recorder) {
        SummaryVals vals = new SummaryVals();
        vals.orderParam = recorder.getMeanOrderParam();
        vals.avgSigmaz = recorder.getMeanAvgZs();
        vals.avgSigmazz = recorder.getMeanAvgZzs();
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
            myWriter.write(SystemParams.getHeader() + ", " + getHeader() + "\n");
            myHeaderOut = true;
        }
        
        // Write the values
        Iterator<SystemParams> iter = myVals.keySet().iterator();
        while(iter.hasNext()) {
            SystemParams params = iter.next();
            SummaryVals vals = myVals.get(params);
            myWriter.write(params.getLine() + ", ");
            myWriter.write(vals.orderParam + ", ");
            myWriter.write(vals.avgSigmaz + ", ");
            myWriter.write(vals.avgSigmazz + "\n");
        }
        
        myWriter.close();
    }
    
    private void writeVals(SystemParams params, SummaryVals vals) {
        if(!myHeaderOut) {
            myWriter.write(SystemParams.getHeader() + ", " + getHeader() + "\n");
            myHeaderOut = true;
        }
        
        myWriter.write(params.getLine() + ", ");
        myWriter.write(vals.orderParam + ", ");
        myWriter.write(vals.avgSigmaz + ", ");
        myWriter.write(vals.avgSigmazz + "\n");
        
        myWriter.flush();
    }
    
    private String getHeader() {
        return "order_param, avg_sigmaz, avg_sigmazz";
    }
}
