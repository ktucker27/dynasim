package utils;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Random;
import java.util.Scanner;

import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;

public class SynchUtils {

    public static void getStartIdx(int[] startIdx, int n) {
        startIdx[0] = 0;
        startIdx[1] = n;
        startIdx[2] = startIdx[1] + n*(n-1);
        startIdx[3] = startIdx[2] + n;
        startIdx[4] = startIdx[3] + n*(n-1)/2;
        startIdx[5] = startIdx[4] + n*(n-1)/2;
    }
    
    public static int getDimension(int n) {
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
    
    public static void initializeHomogeneous(DynaComplex[] z0, double th0, int n) {
        double sth = Math.sin(Math.PI - th0);
        double cth = Math.cos(Math.PI - th0);
        
        int[] startIdx = new int[6];
        getStartIdx(startIdx, n);
        
        for(int i = 0; i < n; ++i) {
            // ps
            z0[i] = new DynaComplex(sth, 0.0);
            
            // zs
            z0[startIdx[2] + i] = new DynaComplex(cth, 0.0);
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
    
    public static void initializeInhomogeneous(DynaComplex[] z0, double th0, int n) {
        int[] startIdx = new int[6];
        getStartIdx(startIdx, n);
        
        for(int i = 0; i < n; ++i) {
            double th = Math.PI/2.0 - th0*Math.cos(7.0*Math.PI*i/6.0);
            
            // ps
            z0[i] = new DynaComplex(Math.sin(th), 0.0);
            
            // zs
            z0[startIdx[2] + i] = new DynaComplex(Math.cos(th), 0.0);
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
    
    public static void detuneLor(double delta, double[] d) {
        int n = d.length;
        
        Random rg = new Random(2);
        double uval;
        for(int i = 0; i < n; ++i) {
            uval = (2.0*rg.nextDouble() - 1.0)*Math.PI*0.5;
            d[i] = delta*Math.tan(uval);
        }
    }
    
    public static DynaComplex compCorr(double[] y0, int n) {
        double sum = 0.0;
        int startIdx = 2*(2*n + n*(n-1));
        
        for(int i = startIdx; i < startIdx + n*(n-1); i += 2) {
            sum += 2.0*y0[i];
        }
        sum *= 1.0/(n*(n-1));
        
        return new DynaComplex(sum);
    }
    
    public static double compSigmazAvg(double[] y0, int n) {
        double sum = 0.0;
        int startIdx = 2*n*n;
        
        for(int i = startIdx; i < startIdx + 2*n; i += 2) {
            sum += y0[i];
        }
        sum *= 1.0/n;

        return sum;
    }
    
    public static DynaComplex compSigmapAvg(double[] y0, int n) {
        DynaComplex sum = new DynaComplex(0.0);
        DynaComplex z = new DynaComplex();
        int startIdx = 0;
        
        for(int i = startIdx; i < startIdx + 2*n; i += 2) {
            z.set(y0[i], y0[i+1]);
            sum.add(z);
        }
        sum.multiply(1.0/n);

        return sum;
    }
    
    public static DynaComplex compSigmazpAvg(double[] y0, int n) {
        int[] startIdxArray = new int[6];
        getStartIdx(startIdxArray, n);
        
        int startIdx = 2*startIdxArray[1];
        int endIdx = 2*startIdxArray[2];
        
        DynaComplex sum = new DynaComplex(0.0);
        DynaComplex z = new DynaComplex();
        
        for(int i = startIdx; i < endIdx; i += 2) {
            z.set(y0[i], y0[i+1]);
            sum.add(z);
        }
        sum.multiply(2.0/(endIdx - startIdx));

        return sum;
    }
    
    public static double compSigmazzAvg(double[] y0, int n) {
        int[] startIdxArray = new int[6];
        getStartIdx(startIdxArray, n);
        
        double sum = 0.0;
        int startIdx = 2*startIdxArray[4];
        int endIdx = 2*startIdxArray[5];
        
        for(int i = startIdx; i < endIdx; i += 2) {
            sum += y0[i];
        }
        sum *= 2.0/(endIdx - startIdx);

        return sum;
    }
    
    public static DynaComplex compSigmappAvg(double[] y0, int n) {
        int[] startIdxArray = new int[6];
        getStartIdx(startIdxArray, n);
        
        int startIdx = 2*startIdxArray[5];
        int endIdx = y0.length;
        
        DynaComplex sum = new DynaComplex(0.0);
        DynaComplex z = new DynaComplex();
        
        for(int i = startIdx; i < endIdx; i += 2) {
            z.set(y0[i], y0[i+1]);
            sum.add(z);
        }
        sum.multiply(2.0/(endIdx - startIdx));

        return sum;
    }

    public static void reduceDim(double[] y1, double[] y2, int n1, int n2) {
        int[] startIdxArray1 = new int[6];
        int[] startIdxArray2 = new int[6];
        getStartIdx(startIdxArray1, n1);
        getStartIdx(startIdxArray2, n2);
        
        // sigma^+/z
        for(int i = 0; i < n2; ++i) {
            y2[2*i] = y1[2*i];
            y2[2*i+1] = y1[2*i+1];
            y2[2*startIdxArray2[2] + 2*i] = y1[2*startIdxArray1[2] + 2*i];
            y2[2*startIdxArray2[2] + 2*i + 1] = y1[2*startIdxArray1[2] + 2*i + 1];
        }
        
        // sigma^+/+/z sigma^+/-/z
        int idx1 = 0;
        int idx2 = 0;
        for(int i = 0; i < n1; ++i) {
            if(i >= n2) break;
            for(int j = i + 1; j < n1; ++j, idx1 += 2) {
                if(j >= n2) continue;
                
                for(int k = 3; k <= 5; ++k) {
                    y2[2*startIdxArray2[k] + idx2] = y1[2*startIdxArray1[k] + idx1];
                    y2[2*startIdxArray2[k] + idx2 + 1] = y1[2*startIdxArray1[k] + idx1 + 1];
                }
                
                idx2 += 2;
            }
        }
        
        // sigma^z sigma^+
        idx1 = 0;
        idx2 = 0;
        for(int i = 0; i < n1; ++i) {
            if(i >= n2) break;
            for(int j = 0; j < n1; ++j, idx1 += 2) {
                if(j >= n2) continue;
                
                if(i == j) {
                    idx1 -= 2;
                    continue;
                }
                
                y2[2*startIdxArray2[1] + idx2] = y1[2*startIdxArray1[1] + idx1];
                y2[2*startIdxArray2[1] + idx2 + 1] = y1[2*startIdxArray1[1] + idx1 + 1];
                
                idx2 += 2;
            }
        }
    }
    
    public static boolean isSteadyState(ArrayList<Double> spinCorr, int numSteps, double tol) {
        boolean isSteadyState = true;
        
        StandardDeviation stdDev = new StandardDeviation();
        double[] zvals = new double[numSteps];
        int idx = 0;
        double sum = 0.0;
        for(int i = spinCorr.size() - numSteps; i < spinCorr.size(); ++i) {
            double corr = spinCorr.get(i);
            
            zvals[idx] = corr;
            sum += corr;
            ++idx;
        }
        
        double mean = sum/(double)numSteps;
        
        if(Math.abs(stdDev.evaluate(zvals)/mean) >= tol) {
            isSteadyState = false;
        }
        
        return isSteadyState;
    }
    
    public static double getW(int n) {
        if(n == 10) {
            return 1.55591;
        }
        
        return 1.515622000024700 + 0.024796919999851*n;
    }
    
    public static double getW_D0(int n) {
        return 1.063413853798558 + 0.029364107193489*n;
    }
    
    public static double getWOpt(int n) {
        return 1.972347641309676 + 0.493905191750637*n;
    }
}
