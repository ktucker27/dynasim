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
    
    public SummaryWriter(String filepath) throws FileNotFoundException {
        myVals = new TreeMap<CumulantParams, SummaryWriter.SummaryVals>();
        myWriter = new PrintWriter(filepath);
    }
    
    public void addVals(SystemEval eval, CumulantParams params, double[] y) {
        SummaryVals vals = new SummaryVals();
        vals.orderParam = eval.getOrderParam(y);
        vals.avgSigmaz = eval.getAvgSigmaz(y);
        myVals.put(params, vals);
    }
    
    public void writeToFile() {
        // Write the header
        myWriter.write(CumulantParams.getHeader() + ", " + getHeader() + "\n");
        
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
    
    private String getHeader() {
        return "order_param, avg_sigmaz";
    }
}
