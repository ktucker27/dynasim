package mcwf;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;

import utils.ExpectedSpinValues;

public class MCWFWriter {
    PrintWriter myWriter;
    
    public MCWFWriter(String filename) throws FileNotFoundException, UnsupportedEncodingException {
        myWriter = new PrintWriter(filename, "UTF-8");
    }
    
    public void write(MCWFAggregator agg) {
        int numTimes = agg.getNumTimes();
        double numTraj = agg.getNumTrajectories();
        
        for(int timeIdx = 0; timeIdx < numTimes; ++timeIdx) {
            ExpectedSpinValues evs = agg.getEvs(timeIdx);
            
            myWriter.print(agg.getTime(timeIdx));
            
            for(int i = 0; i < 3; ++i) {
                myWriter.print(", " + evs.getEs(i)/numTraj);
            }
            
            for(int i = 0; i < 3; ++i) {
                for(int j = 0; j < 3; ++j) {
                    myWriter.print(", " + evs.getEss(i, j).getReal()/numTraj);
                    myWriter.print(", " + evs.getEss(i, j).getImaginary()/numTraj);
                }
            }
            myWriter.print("\n");
        }
        
        myWriter.close();
    }
}
