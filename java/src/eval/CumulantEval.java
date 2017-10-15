package eval;

import java.util.Random;

import ode.DynaComplexODEAdapter;
import utils.DynaComplex;


public class CumulantEval implements SystemEval {

    int n;
    int[] startIdx;
    DynaComplex[] z;
    
    public CumulantEval(int n) {
        this.n = n;

        startIdx = new int[6];
        startIdx[0] = 0;
        startIdx[1] = n;
        startIdx[2] = startIdx[1] + n*(n-1);
        startIdx[3] = startIdx[2] + n;
        startIdx[4] = startIdx[3] + n*(n-1)/2;
        startIdx[5] = startIdx[4] + n*(n-1)/2;
        
        z = new DynaComplex[getDimension()];
        for(int i = 0; i < z.length; ++i) {
            z[i] = new DynaComplex();
        }
    }
    
    public int[] getStartIdx() {
        return startIdx;
    }
    
    @Override
    public int getN() {
        return n;
    }
    
    @Override
    public int getRealDimension() {
        return 2*getDimension();
    }

    @Override
    public int getDimension() {
        return n*(n-1) + 3*n*(n-1)/2 + 2*n;
    }

    @Override
    public double getOrderParam(double[] y) {
        double sum = 0.0;
        int startIdx = 2*(2*n + n*(n-1));
        
        for(int i = startIdx; i < startIdx + n*(n-1); i += 2) {
            sum += 2.0*y[i];
        }
        sum *= 1.0/(n*(n-1));
        
        return sum;
    }
    
    @Override
    public double getAvgSigmaz(double[] y) {
        double sum = 0.0;
        int beginIdx = 2*startIdx[2];
        
        for(int i = beginIdx; i < beginIdx + 2*n; i += 2) {
            sum += y[i];
        }
        sum *= 1.0/n;

        return sum;
    }

    @Override
    public void getBlochVectors(double[] y, double[] xs, double[] ys, double[] zs) {
        // TODO
        throw new UnsupportedOperationException("Bloch vector retrieval not yet implemented for CumulantEval");
    }
    
    @Override
    public void initSpinUpX(double[] y0) {
        initialize(y0, 0.5*Math.PI, 0.0, InitAngleType.CONST, InitAngleType.CONST);
    }
    
    @Override
    public void initialize(double[] y0, double zenith, double phase, InitAngleType zenithType, InitAngleType phaseType) {
        double zi = 0.0;
        double pi = 0.0;
        double szi = 0.0;
        double czi = 0.0;
        
        double eps = 1.0e-3;
        
//        Scanner inputStream = new Scanner(new File("/Users/kristophertucker/Google Drive/Research/Synch/cumulant_all/phase.txt"));
//        int idx = 0;
//        while(inputStream.hasNext()) {
//            phase[idx] = Double.parseDouble(inputStream.next());
//            ++idx;
//            if(idx >= n) break;
//        }
        
        Random rg = new Random(1);
        for(int i = 0; i < n; ++i) {
            switch(zenithType) {
            case EQUAL_SPACING:
                zi = eps + i*(Math.PI - 2*eps)/(double)(n-1);
                break;
            case RANDOM:
                zi = rg.nextDouble()*Math.PI;
                break;
            case CONST:
                zi = zenith;
                break;
            }
            
            szi = Math.sin(zi);
            czi = Math.cos(zi);
            
            switch(phaseType) {
            case EQUAL_SPACING:
                pi = i*2*Math.PI/(double)n - Math.PI;
                break;
            case RANDOM:
                pi = rg.nextDouble()*2.0*Math.PI - Math.PI;
                break;
            case CONST:
                pi = phase;
                break;
            }
            
            // ps
            z[i] = new DynaComplex(0.5*szi*Math.cos(pi),
                                    0.5*szi*Math.sin(pi));
            
            // zs
            z[startIdx[2] + i] = new DynaComplex(czi, 0.0);
        }
        
        int idx = 0;
        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < n; ++j) {
                if(i == j) {
                    continue;
                }
                
                // zps
                z[startIdx[1] + idx] = new DynaComplex(z[startIdx[2] + i]);
                z[startIdx[1] + idx].multiply(z[j]);
                ++idx;
            }
        }
        
        idx = 0;
        DynaComplex temp = new DynaComplex();
        for(int i = 0; i < n; ++i) {
            for(int j = i + 1; j < n; ++j) {
                // pms
                z[startIdx[3] + idx] = new DynaComplex(z[i]);
                z[startIdx[3] + idx].multiply(temp.set(z[j]).conjugate());
                
                // zzs
                z[startIdx[4] + idx] = new DynaComplex(z[startIdx[2] + i]);
                z[startIdx[4] + idx].multiply(z[startIdx[2] + j]);
                
                // pps
                z[startIdx[5] + idx] = new DynaComplex(z[i]);
                z[startIdx[5] + idx].multiply(z[j]);
                
                ++idx;
            }
        }
        
        DynaComplexODEAdapter.toReal(z, y0);
    }
}
