package eval;

import java.util.Random;

import ode.DynaComplexODEAdapter;
import utils.DynaComplex;


public class CumulantEval implements SystemEval {

    int n;
    int[] startIdx;
    DynaComplex[] z;
    
    DynaComplex ct1, ct2, ct3;
    
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
        
        ct1 = new DynaComplex();
        ct2 = new DynaComplex();
        ct3 = new DynaComplex();
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
    public double getAvgSigmazz(double[] y) {
        double sum = 0.0;
        int beginIdx = 2*startIdx[4];
        int endIdx = 2*startIdx[5];
        
        for(int i = beginIdx; i < endIdx; i += 2) {
            sum += y[i];
        }
        sum *= 1.0/(startIdx[5] - startIdx[4]);

        return sum;
    }

    @Override
    public void getBlochVectors(double[] y, double[] xs, double[] ys, double[] zs) {
        for(int i = 0; i < n; ++i) {
            xs[i] = 2*y[2*i];
            ys[i] = 2*y[2*i+1];
            zs[i] = y[2*startIdx[2] + 2*i];
        }
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
    
    private int getTriIdx(int i, int j) {
        if(i > j) {
            int t = j;
            j = i;
            i = t;
        }
        
        return n*i - i*(i+1)/2 + j - i - 1;
    }
    
    private int getRecIdx(int i, int j) {
        if(j < i) {
            return i*(n-1) + j;
        }
        
        return i*(n-1) + j - 1;
    }
    
    // In the following, for the superscript indices:
    // 0 = +, 1 = -, 2 = z
    
    private enum CumulantOp {
        PLUS,
        MINUS,
        Z,
        ZERO,
        IDENTITY
    }
    
    private class ReducedSingle {
        public CumulantOp op;
        public double m;
        public double b;
        
        public ReducedSingle() {
            op = CumulantOp.IDENTITY;
            m = 1.0;
            b = 0.0;
        }
    }
    
    private ReducedSingle getSingleOp(int al, int bt) {
        ReducedSingle single = new ReducedSingle();
        if(al == bt) {
            if(al == 0 || al == 1) {
                single.op = CumulantOp.ZERO;
                single.m = 0.0;
            } else {
                single.op = CumulantOp.IDENTITY;
                single.m = 1.0;
            }
        } else {
            switch(al) {
                case 0:
                    if(bt == 1) {
                        single.op = CumulantOp.Z;
                        single.m = 0.5;
                        single.b = 0.5;
                    } else {
                        single.op = CumulantOp.PLUS;
                        single.m = -1.0;
                        single.b = 0.0;
                    }
                    break;
                case 1:
                    if(bt == 0) {
                        single.op = CumulantOp.Z;
                        single.m = -0.5;
                        single.b = 0.5;
                    } else {
                        single.op = CumulantOp.MINUS;
                        single.m = 1.0;
                        single.b = 0.0;
                    }
                    break;
                case 2:
                    if(bt == 0) {
                        single.op = CumulantOp.PLUS;
                        single.m = 1.0;
                        single.b = 0.0;
                    } else {
                        single.op = CumulantOp.MINUS;
                        single.m = -1.0;
                        single.b = 0.0;
                    }
                    break;
                default:
                    throw new UnsupportedOperationException("Invalid alpha passed to getSingleOp");
            }
        }
        
        return single;
    }
    
    /**
     * @param al superscript
     * @param a particle index
     * @param z value array
     * @return <sigma_a^al> taken from z
     */
    public DynaComplex getSingle(int al, int a, DynaComplex[] z, DynaComplex ans) {
        ans.set(z[a + startIdx[2]*(al/2)]);
        if(al == 1) ans.conjugate();
        return ans;
    }
    
    /**
     * @param al first superscript
     * @param bt second superscript
     * @param a first particle index
     * @param b second particle index
     * @param z value array
     * @return <sigma_a^al sigma_b^bt> taken from z
     */
    public DynaComplex getDouble(int al, int bt, int a, int b, DynaComplex[] z, DynaComplex ans) {
        // Handle same particle numbers first
        if(a == b) {
            ReducedSingle single = getSingleOp(al, bt);
            
            if(single.op == CumulantOp.ZERO) {
                ans.set(0, 0);
            } else if(single.op == CumulantOp.IDENTITY){
                ans.set(1, 0);
            } else {
                getSingle(single.op.ordinal(), a, z, ans);
                ans.multiply(single.m).add(single.b);
            }
            
            return ans;
        }
        
        switch(al) {
        case 0:
            switch(bt) {
            case 0:
                return ans.set(z[startIdx[5] + getTriIdx(a,b)]);
            case 1:
                ans.set(z[startIdx[3] + getTriIdx(a,b)]);
                if(a > b) return ans.conjugate();
                return ans;
            case 2:
                return ans.set(z[startIdx[1] + getRecIdx(b,a)]);
            }
        case 1:
            switch(bt) {
            case 0:
                ans.set(z[startIdx[3] + getTriIdx(a,b)]);
                if(a < b) return ans.conjugate();
                return ans;
            case 1:
                return ans.set(z[startIdx[5] + getTriIdx(a,b)]).conjugate();
            case 2:
                return ans.set(z[startIdx[1] + getRecIdx(b,a)]).conjugate();
            }
        case 2:
            switch(bt) {
            case 0:
                return ans.set(z[startIdx[1] + getRecIdx(a,b)]);
            case 1:
                return ans.set(z[startIdx[1] + getRecIdx(a,b)]).conjugate();
            case 2:
                return ans.set(z[startIdx[4] + getTriIdx(a,b)]);
            }
        }
        
        // This should never happen
        throw new UnsupportedOperationException("Invalid superscripts passed to cumulantDouble: " + al + " " + bt);
    }
    
    /**
     * Returns <sigma_a^al sigma_b^bt sigma_c^gm> approximated using the cumulant expansion
     * 
     * @param al first superscript
     * @param bt second superscript
     * @param gm third superscript
     * @param a first particle index
     * @param b second particle index
     * @param c third particle index
     * @param z value array
     * @return <sigma_a^al sigma_b^bt sigma_c^gm> taken from z using the cumulant expansion
     */
    public DynaComplex getTriple(int al, int bt, int gm, int a, int b, int c, DynaComplex[] z, DynaComplex ans) {
        if(a == c) {
            ReducedSingle single = getSingleOp(al, gm);
            if(single.op == CumulantOp.ZERO) {
                ans.set(0, 0);
            } else if(single.op == CumulantOp.IDENTITY) {
                getSingle(bt, b, z, ans);
            } else {
                getDouble(single.op.ordinal(), bt, a, b, z, ans);
                ans.multiply(single.m).add(getSingle(bt, b, z, ct1).multiply(single.b));
            }
        } else if(b == c) {
            ReducedSingle single = getSingleOp(bt, gm);
            if(single.op == CumulantOp.ZERO) {
                ans.set(0, 0);
            } else if(single.op == CumulantOp.IDENTITY) {
                getSingle(al, a, z, ans);
            } else {
                getDouble(al, single.op.ordinal(), a, b, z, ans);
                ans.multiply(single.m).add(getSingle(al, a, z, ct1).multiply(single.b));
            }
        } else {
            getSingle(gm, c, z, ans).multiply(getDouble(al, bt, a, b, z, ct1));
            ans.add(getSingle(al, a, z, ct1).multiply(getDouble(bt, gm, b, c, z, ct2)));
            ans.add(getSingle(bt, b, z, ct1).multiply(getDouble(al, gm, a, c, z, ct2)));
            ans.subtract(getSingle(al, a, z, ct1).multiply(getSingle(bt, b, z, ct2)).multiply(getSingle(gm, c, z, ct3)).multiply(2.0));
        }
        
        return ans;
    }
}
