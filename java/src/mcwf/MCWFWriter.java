package mcwf;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.Iterator;
import java.util.Map;

import utils.DynaComplex;
import utils.ExpectedSpinValues;
import utils.HusimiDist;

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
    
    public void writeAllStates(QuantumTrajectory[] trajectories, String filename) throws FileNotFoundException, UnsupportedEncodingException {
        PrintWriter writer = new PrintWriter(filename, "UTF-8");
        
        for(int trajIdx = 0; trajIdx < trajectories.length; ++trajIdx) {
            writer.print(trajectories[trajIdx].getJ());
            
            DynaComplex[] state = trajectories[trajIdx].getState();
            for(int valIdx = 0; valIdx < state.length; ++valIdx) {
                writer.print(", " + state[valIdx].getReal() + ", " + state[valIdx].getImaginary());
            }
            
            writer.print("\n");
        }
        
        writer.close();
    }
    
    public void writeHusimi(MCWFAggregator agg, String outdir) throws FileNotFoundException, UnsupportedEncodingException {
        Map<Double, HusimiDist> husimis = agg.getHusimis();
        Iterator<Map.Entry<Double, HusimiDist>> iter = husimis.entrySet().iterator();
        while(iter.hasNext()) {
            Map.Entry<Double, HusimiDist> entry = iter.next();
            String timeStr = String.format("%.2f", entry.getKey()).replace('.', 'p');
            String filename = outdir + "/husimi_" + timeStr + ".txt";
            PrintWriter writer = new PrintWriter(filename, "UTF-8");
            writer.print(agg.getNumTrajectories() + "\n");
            writeHusimi(entry.getValue(), writer);
            writer.close();
        }
    }
    
    public void writeHusimi(HusimiDist dist, PrintWriter writer) {
        if(dist.isZero()) {
            writer.print("0\n");
            return;
        }
        
        final double[][] vals = dist.getVals();
        for(int i = 0; i < dist.getNumTheta(); ++i) {
            for(int j = 0; j < dist.getNumPhi(); ++j) {
                if(j != 0) writer.print(", ");
                writer.print(vals[i][j]);
            }
            writer.print("\n");
        }
    }
}
