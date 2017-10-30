package utils;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Random;
import java.util.Scanner;
import java.util.concurrent.atomic.AtomicInteger;

import org.apache.commons.math3.ode.nonstiff.AdamsMoultonIntegrator;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;

import coupling.DynaConstCoupling;
import handlers.TwoTimeHandler;
import integrator.IntegratorRequest;
import ode.CorrelationODEs;
import ode.CumulantParams;
import ode.DynaComplexODEAdapter;
import ode.FOCorrelationODEs;

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
    
    public static void initializeConst(DynaComplex[] z0, double tip, int n) {
        int[] startIdx = new int[6];
        getStartIdx(startIdx, n);
        
        for(int i = 0; i < n; ++i) {
            // ps
            z0[i] = new DynaComplex(0.5, 0.0);
            
            // zs
            z0[startIdx[2] + i] = new DynaComplex(0.0, 0.0);
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
        inputStream.close();
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
    
    public static void detuneDiscrete(double delta, double[] d) {
        int n = d.length;
        
        for(int i = 0; i < n; ++i) {
            if(i % 2 == 0) {
                d[i] = delta;
            } else {
                d[i] = -1.0*delta;
            }
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

    public static void changeDim(double[] y1, double[] y2, int n1, int n2) {
        int[] startIdxArray1 = new int[6];
        int[] startIdxArray2 = new int[6];
        getStartIdx(startIdxArray1, n1);
        getStartIdx(startIdxArray2, n2);
        
        int nmin, nmax;
        if(n1 < n2) {
            nmin = n1;
            nmax = n2;
        } else {
            nmin = n2;
            nmax = n1;
        }
        
        // sigma^+/z
        for(int i = 0; i < nmin; ++i) {
            y2[2*i] = y1[2*i];
            y2[2*i+1] = y1[2*i+1];
            y2[2*startIdxArray2[2] + 2*i] = y1[2*startIdxArray1[2] + 2*i];
            y2[2*startIdxArray2[2] + 2*i + 1] = y1[2*startIdxArray1[2] + 2*i + 1];
        }
        
        // sigma^+/+/z sigma^+/-/z
        AtomicInteger idx1 = new AtomicInteger(0);
        AtomicInteger idx2 = new AtomicInteger(0);
        AtomicInteger largerIdx;
        AtomicInteger smallerIdx;
        if(n1 < n2) {
            largerIdx = idx2;
            smallerIdx = idx1;
        } else {
            largerIdx = idx1;
            smallerIdx = idx2;
        }
        
        for(int i = 0; i < nmax; ++i) {
            if(i >= nmin) break;
            for(int j = i + 1; j < nmax; ++j, largerIdx.getAndAdd(2)) {
                if(j >= nmin) continue;
                
                for(int k = 3; k <= 5; ++k) {
                    y2[2*startIdxArray2[k] + idx2.get()] = y1[2*startIdxArray1[k] + idx1.get()];
                    y2[2*startIdxArray2[k] + idx2.get() + 1] = y1[2*startIdxArray1[k] + idx1.get() + 1];
                }
                
                smallerIdx.getAndAdd(2);
            }
        }
        
        // sigma^z sigma^+
        idx1.set(0);
        idx2.set(0);
        for(int i = 0; i < nmax; ++i) {
            if(i >= nmin) break;
            for(int j = 0; j < nmax; ++j, largerIdx.getAndAdd(2)) {
                if(j >= nmin) continue;
                
                if(i == j) {
                    largerIdx.getAndAdd(-2);
                    continue;
                }
                
                y2[2*startIdxArray2[1] + idx2.get()] = y1[2*startIdxArray1[1] + idx1.get()];
                y2[2*startIdxArray2[1] + idx2.get() + 1] = y1[2*startIdxArray1[1] + idx1.get() + 1];
                
                smallerIdx.getAndAdd(2);
            }
        }
        
        // If we are increasing the dimension, we need to go through the new rows and
        // columns and make sure the values are consistent
        if(n2 > n1) {
            DynaComplex z = new DynaComplex();
            DynaComplex t1 = new DynaComplex();
            
            Random rg = new Random(0);
            for(int i = n1; i < n2; ++i) {
                double phase = rg.nextDouble()*2.0*Math.PI;
                
                // ps
                y2[2*i] = 0.5*Math.cos(phase);
                y2[2*i + 1] = 0.5*Math.sin(phase);
                
                // zs
                y2[2*startIdxArray2[2] + 2*i] = 0.0;
                y2[2*startIdxArray2[2] + 2*i + 1] = 0.0;
            }
            
            int idx = 0;
            for(int i = 0; i < n2; ++i) {
                for(int j = 0; j < n2; ++j) {
                    if(i == j) {
                        continue;
                    }
                    
                    // zps
                    if(i >= n1 || j >= n1) {
                        z.set(y2[2*startIdxArray2[2] + 2*i], y2[2*startIdxArray2[2] + 2*i + 1]);
                        z.multiply(t1.set(y2[2*j], y2[2*j + 1]));
                        y2[2*startIdxArray2[1] + 2*idx] = z.getReal();
                        y2[2*startIdxArray2[1] + 2*idx + 1] = z.getImaginary();
                    }
                    ++idx;
                }
            }
            
            idx = 0;
            for(int i = 0; i < n2; ++i) {
                for(int j = i + 1; j < n2; ++j) {
                    if(i >= n1 || j >= n1) {
                        // pms
                        z.set(y2[2*i], y2[2*i + 1]);
                        z.multiply(t1.set(y2[2*j], y2[2*j + 1]).conjugate());
                        y2[2*startIdxArray2[3] + 2*idx] = z.getReal();
                        y2[2*startIdxArray2[3] + 2*idx + 1] = z.getImaginary();

                        // zzs
                        z.set(y2[2*startIdxArray2[2] + 2*i], y2[2*startIdxArray2[2] + 2*i + 1]);
                        z.multiply(t1.set(y2[2*startIdxArray2[2] + 2*j], y2[2*startIdxArray2[2] + 2*j + 1]));
                        y2[2*startIdxArray2[4] + 2*idx] = z.getReal();
                        y2[2*startIdxArray2[4] + 2*idx + 1] = z.getImaginary();

                        // pps
                        z.set(y2[2*i], y2[2*i + 1]);
                        z.multiply(t1.set(y2[2*j], y2[2*j + 1]));
                        y2[2*startIdxArray2[5] + 2*idx] = z.getReal();
                        y2[2*startIdxArray2[5] + 2*idx + 1] = z.getImaginary();
                    }
                    ++idx;
                }
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
    
    private static int getTriIdx(int i, int j, int n) {
        if(i > j) {
            int t = j;
            j = i;
            i = t;
        }
        
        return n*i - i*(i+1)/2 + j - i - 1;
    }
    
    private static int getRecIdx(int i, int j, int n) {
        if(j < i) {
            return i*(n-1) + j;
        }
        
        return i*(n-1) + j - 1;
    }
    
    public static IntegratorRequest getCorrRequest(CumulantParams params, double[] y, String filename) throws FileNotFoundException, UnsupportedEncodingException {
        int n = params.getN();
        
        int dim = SynchUtils.getDimension(n);
        int[] startIdx = new int[6];
        SynchUtils.getStartIdx(startIdx, n);
        DynaComplex[] z = new DynaComplex[dim];
        for(int i = 0; i < dim; ++i) {
            z[i] = new DynaComplex(0, 0);
        }
        DynaComplexODEAdapter.toComplex(y, z);
        
        DynaComplex[] szs = new DynaComplex[n];
        int idx2 = 0;
        for(int i = startIdx[2]; i < startIdx[2] + n; ++i) {
            szs[idx2] = new DynaComplex(z[i].getReal(), z[i].getImaginary());
            ++idx2;
        }
        
        DynaComplex[] z02 = new DynaComplex[n*n];
        int idx1 = startIdx[2];
        idx2 = startIdx[3];
        for(int i = 0; i < n; ++i) {
            z02[i*n + i] = new DynaComplex(0.5 + 0.5*z[idx1].getReal(), 0);
            ++idx1;
            for(int j = i+1; j < n; ++j) {
                z02[i*n + j] = new DynaComplex(z[idx2].getReal(), z[idx2].getImaginary());
                z02[j*n + i] = new DynaComplex(z[idx2].getReal(), z[idx2].getImaginary());
                z02[j*n + i].conjugate();
                
                ++idx2;
            }
        }

        DynaConstCoupling coupling = new DynaConstCoupling(params.getAlpha().getReal(), params.getAlpha().getImaginary());
        CorrelationODEs c_corr_odes = new CorrelationODEs(n, params.getGamma(), params.getW(), coupling, params.getD(), szs);
        DynaComplexODEAdapter odes = new DynaComplexODEAdapter(c_corr_odes);
        
        double[] y02 = new double[2*n*n];
        DynaComplexODEAdapter.toReal(z02, y02);
        
        FullODESolution soln = new FullODESolution();
        IntegratorRequest request = new IntegratorRequest(odes, 0, y02, 5, soln);
        
        if(!filename.isEmpty()) {
//          int[] out_col = {0,1,2,3,2*(n+1),2*(n+1)+1};
//          WriteHandler writeHandler = new WriteHandler(filename, out_col);
            TwoTimeHandler writeHandler = new TwoTimeHandler(filename, n, false);
            request.addStepHandler(writeHandler);
        }
        
        return request;
    }
    
    public static IntegratorRequest getFOCorrRequest(CumulantParams params, double[] y, String filename) throws FileNotFoundException, UnsupportedEncodingException
    {
        int n = params.getN();
        
        int dim = SynchUtils.getDimension(n);
        int[] startIdx = new int[6];
        SynchUtils.getStartIdx(startIdx, n);
        DynaComplex[] z = new DynaComplex[dim];
        for(int i = 0; i < dim; ++i) {
            z[i] = new DynaComplex(0, 0);
        }
        DynaComplexODEAdapter.toComplex(y, z);
        
        DynaComplex[] szs = new DynaComplex[n];
        int idx2 = 0;
        for(int i = startIdx[2]; i < startIdx[2] + n; ++i) {
            szs[idx2] = new DynaComplex(z[i].getReal(), z[i].getImaginary());
            ++idx2;
        }
        
        DynaComplex[] z02 = new DynaComplex[n*n + n*n*(n-1)];
        int idx1 = startIdx[2];
        idx2 = startIdx[3];
        for(int i = 0; i < n; ++i) {
            z02[i*n + i] = new DynaComplex(0.5 + 0.5*z[idx1].getReal(), 0);
            ++idx1;
            for(int j = i+1; j < n; ++j) {
                z02[i*n + j] = new DynaComplex(z[idx2].getReal(), z[idx2].getImaginary());
                z02[j*n + i] = new DynaComplex(z[idx2].getReal(), z[idx2].getImaginary());
                z02[j*n + i].conjugate();
                
                ++idx2;
            }
        }
        
        int idx = n*n;
        for(int c = 0; c < n; ++c) {
            for(int a = 0; a < n; ++a) {
                for(int b = 0; b < n; ++b) {
                    if(a == b) continue;
                    
                    z02[idx] = new DynaComplex();
                    if(c == a) {
                        z02[idx].set(z[startIdx[3] + getTriIdx(a,b,n)]);
                        z02[idx].conjugate();
                    } else if(c == b) {
                        z02[idx].set(0.5*(szs[a].getReal() + z[startIdx[4] + getRecIdx(a,b,n)].getReal()), 0);
                    } else {
                        z02[idx].set(z[startIdx[3] + getTriIdx(b,c,n)]).multiply(szs[a]);
                    }
                    
                    ++idx;
                }
            }
        }
        
//        for(int i = 0; i < z02.length; ++i) {
//            System.out.println(z02[i]);
//        }

        DynaConstCoupling coupling = new DynaConstCoupling(params.getAlpha().getReal(), params.getAlpha().getImaginary());
        FOCorrelationODEs c_corr_odes = new FOCorrelationODEs(n, params.getGamma(), params.getW(), coupling, params.getD(), szs);
        DynaComplexODEAdapter odes = new DynaComplexODEAdapter(c_corr_odes);
        
        double[] y02 = new double[2*z02.length];
        DynaComplexODEAdapter.toReal(z02, y02);
        
        FullODESolution soln = new FullODESolution();
        IntegratorRequest request = new IntegratorRequest(odes, 0, y02, 5, soln);
        
        if(!filename.isEmpty()) {
//          int[] out_col = {0,1,2,3,2*(n+1),2*(n+1)+1};
//          WriteHandler writeHandler = new WriteHandler(filename, out_col);
            TwoTimeHandler writeHandler = new TwoTimeHandler(filename, n, false);
            request.addStepHandler(writeHandler);
        }
        
        return request;
    }
    
    public static double[] compCorr(CumulantParams params, double[] y, String filename, boolean fo) throws FileNotFoundException, UnsupportedEncodingException
    {
        IntegratorRequest request;
        if(fo) {
            request = getFOCorrRequest(params, y, filename);
        } else {
            request = getCorrRequest(params, y, filename);
        }
        
        AdamsMoultonIntegrator integrator = new AdamsMoultonIntegrator(2, 1.0e-18, .0001, 1.0e-3, 1.0e-2);
        for(int i = 0; i < request.numStepHandlers(); ++i) {
            integrator.addStepHandler(request.getStepHandler(i));
        }

        for(int i = 0; i < request.numEventHandlers(); ++i) {
            integrator.addEventHandler(request.getEventHandler(i), Double.POSITIVE_INFINITY, 1.0e-12, 100);
        }
        
        double[] y2 = new double[request.getY0().length];
        
        double startTime = System.nanoTime();
        
        integrator.integrate(request.getOdes(), request.getT0(), request.getY0(), request.getTF(), y2);
        request.getSoln().setSolution(y2);
        
        double endTime = System.nanoTime();
        
        System.out.println("Correlation time: " + (endTime - startTime)/1.0e9 + " seconds");
        
        return y2;
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
    
    public static int pow(int base, int exp) {
        int ans = 1;
        for(int i = 0; i < exp; ++i) {
            ans *= base;
        }
        
        return ans;
    }
}
