package utils;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Random;
import java.util.Scanner;

public class SynchUtils {

    public static void getStartIdx(int[] startIdx, int n) {
        startIdx[0] = 0;
        startIdx[1] = n;
        startIdx[2] = startIdx[1] + n*(n-1);
        startIdx[3] = startIdx[2] + n;
        startIdx[4] = startIdx[3] + n*(n-1)/2;
        startIdx[5] = startIdx[4] + n*(n-1)/2;
    }
    
    public static int getDim(int n) {
        return n*(n-1) + 3*n*(n-1)/2 + 2*n;
    }
    
    public static void initialize(DynaComplex[] z0, double tip, int n) {
        double stip2 = Math.sin(tip/2.0);
        
        int[] startIdx = new int[6];
        getStartIdx(startIdx, n);
        
        double[] phase = new double[n];
        
//        Scanner inputStream = new Scanner(new File("/Users/kristophertucker/Google Drive/Research/Synch/cumulant_all/phase.txt"));
//        int idx = 0;
//        while(inputStream.hasNext()) {
//            phase[idx] = Double.parseDouble(inputStream.next());
//            ++idx;
//            if(idx >= n) break;
//        }
        
        Random rg = new Random(1);
        for(int i = 0; i < n; ++i) {
            phase[i] = rg.nextDouble()*2.0*Math.PI;
            
            // ps
            z0[i] = new DynaComplex(stip2*stip2*Math.cos(phase[i]),
                                    stip2*stip2*Math.sin(phase[i]));
            
            // zs
            z0[startIdx[2] + i] = new DynaComplex(Math.cos(tip), 0.0);
        }
        
        int idx = 0;
        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < n; ++j) {
                if(i == j) {
                    continue;
                }
                
                // zps
                z0[startIdx[1] + idx] = new DynaComplex(z0[startIdx[2] + i]);
                z0[startIdx[1] + idx].multiply(z0[j]);
                ++idx;
            }
        }
        
        idx = 0;
        DynaComplex temp = new DynaComplex();
        for(int i = 0; i < n; ++i) {
            for(int j = i + 1; j < n; ++j) {
                // pms
                z0[startIdx[3] + idx] = new DynaComplex(z0[i]);
                z0[startIdx[3] + idx].multiply(temp.set(z0[j]).conjugate());
                
                // zzs
                z0[startIdx[4] + idx] = new DynaComplex(z0[startIdx[2] + i]);
                z0[startIdx[4] + idx].multiply(z0[startIdx[2] + j]);
                
                // pps
                z0[startIdx[5] + idx] = new DynaComplex(z0[i]);
                z0[startIdx[5] + idx].multiply(z0[j]);
                
                ++idx;
            }
        }
    }
    
    public static void detuneFile(String filename, double[] d) throws FileNotFoundException {
        int n = d.length;
        Scanner inputStream = new Scanner(new File(filename));
        int idx = 0;
        while(inputStream.hasNext()) {
            d[idx] = Double.parseDouble(inputStream.next());
            ++idx;
            if(idx >= n) break;
        }
    }
    
    public static void detuneGauss(double delta, double[] d) {
        int n = d.length;
        
        Random rg = new Random(2);
        for(int i = 0; i < n; ++i) {
            d[i] = delta*rg.nextGaussian();
        }
    }
}
