package mcwf;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;

import utils.ExpectedSpinValues;

public class MCWFWriter {
    
    public MCWFWriter() {
    }
    
    public void write(MCWFAggregator agg, String filename) throws FileNotFoundException, UnsupportedEncodingException {
        PrintWriter writer = new PrintWriter(filename, "UTF-8");

        int numTimes = agg.getNumTimes();
        
        // Print the header
        writer.print(agg.getNumTrajectories() + "\n");
        
        for(int timeIdx = 0; timeIdx < numTimes; ++timeIdx) {
            ExpectedSpinValues evs = agg.getSumEvs(timeIdx);
            
            writer.print(agg.getTime(timeIdx));
            
            // Write sums
            for(int i = 0; i < 3; ++i) {
                writer.print(", " + evs.getEs(i));
            }
            
            for(int i = 0; i < 3; ++i) {
                for(int j = 0; j < 3; ++j) {
                    writer.print(", " + evs.getEss(i, j).getReal());
                    writer.print(", " + evs.getEss(i, j).getImaginary());
                }
            }
            
            // Write sum squares
            evs = agg.getSumSqEvs(timeIdx);
            for(int i = 0; i < 3; ++i) {
                writer.print(", " + evs.getEs(i));
            }
            
            for(int i = 0; i < 3; ++i) {
                for(int j = 0; j < 3; ++j) {
                    writer.print(", " + evs.getEss(i, j).getReal());
                    writer.print(", " + evs.getEss(i, j).getImaginary());
                }
            }
            
            writer.print("\n");
        }
        
        writer.close();
    }
    
    public void writeJumps(MCWFAggregator agg, String filename) throws FileNotFoundException, UnsupportedEncodingException {
        PrintWriter writer = new PrintWriter(filename, "UTF-8");

        int numTimes = agg.getNumTimes();
        
        // Print the header
        writer.print(agg.getNumTrajectories() + "\n");
        
        for(int timeIdx = 0; timeIdx < numTimes; ++timeIdx) {
            int[] jumps = agg.getJumps(timeIdx);
            
            writer.print(agg.getTime(timeIdx));
            
            // Write jumps
            for(int i = 0; i < jumps.length; ++i) {
                writer.print(", " + jumps[i]);
            }
            
            writer.print("\n");
        }
        
        writer.close();
    }
    
    public void writeState(MCWFAggregator agg, String filename) throws FileNotFoundException, UnsupportedEncodingException {
        PrintWriter writer = new PrintWriter(filename, "UTF-8");

        int numTimes = agg.getNumTimes();
        
        // Print the header
        writer.print(agg.getNumTrajectories() + "\n");
        
        for(int timeIdx = 0; timeIdx < numTimes; ++timeIdx) {
            writer.print(agg.getTime(timeIdx));
            
            writer.print(", " + agg.getSumJs(timeIdx));
            
            writer.print("\n");
        }
        
        writer.close();
    }
}
