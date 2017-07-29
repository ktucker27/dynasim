package handlers;

import java.io.ByteArrayOutputStream;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;

import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.sampling.StepHandler;
import org.apache.commons.math3.ode.sampling.StepInterpolator;

public class WriteHandler implements StepHandler {
    PrintWriter myWriter;
    ByteArrayOutputStream myBuffer;
    String myFilename;
    int[] myWriteIdx;
    
    public WriteHandler(String filename, int[] writeIdx, boolean liveStream) throws FileNotFoundException, UnsupportedEncodingException {
        myWriteIdx = writeIdx;
        myFilename = filename;
        
        if(liveStream) {
            myBuffer = null;
            myWriter = new PrintWriter(filename, "UTF-8");
        } else {
            myBuffer = new ByteArrayOutputStream();
            myWriter = new PrintWriter(myBuffer);
        }
    }
    
    @Override
    public void init(double t0, double[] y0, double t) {
        myWriter.print("0");
        for(int i = 0; i < myWriteIdx.length; ++i) {
            myWriter.print(", " + y0[myWriteIdx[i]]);
        }
        myWriter.print("\n");
    }
    
    @Override
    public void handleStep(StepInterpolator interpolator, boolean isLast)
            throws MaxCountExceededException {
        double[] y = interpolator.getInterpolatedState();
        
        myWriter.print(interpolator.getInterpolatedTime());
        for(int i = 0; i < myWriteIdx.length; ++i) {
            myWriter.print(", " + y[myWriteIdx[i]]);
        }
        myWriter.print("\n");
        
        if(!isLast) {
            myWriter.flush();
        } else {
            myWriter.close();
            if(myBuffer != null) {
                try {
                    PrintWriter fileWriter = new PrintWriter(myFilename, "UTF-8");
                    fileWriter.print(myBuffer.toString());
                    fileWriter.close();
                } catch(FileNotFoundException ex) {
                    System.err.println("Could not open output file " + myFilename);
                    System.err.println(ex.getMessage());
                } catch(UnsupportedEncodingException ex) {
                    System.err.println(ex.getMessage());
                }
            }
        }
    }
}
