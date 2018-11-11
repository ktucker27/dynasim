package eval;

import utils.DynaComplex;

public class SymmEval implements SystemEval {
    int n;
    int [][][] myIdxMap;
    
    public SymmEval(int n) {
        this.n = n;
        myIdxMap = getIdxMap();
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
        return (n+3)*(n+2)*(n+1)/6;
    }

    @Override
    public double getOrderParam(double[] y) {
        // TODO
        return 0;
    }

    @Override
    public double getAvgSigmaz(double[] y) {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public double getAvgSigmazz(double[] y) {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public void getBlochVectors(double[] y, double[] xs, double[] ys, double[] zs) {
        // TODO Auto-generated method stub

    }
    
    @Override
    public void getFirstOrderCollectiveEvs(double[] y, double[] es) {
        int idx001 = myIdxMap[0][0][1];
        int idx010 = myIdxMap[0][1][0];
        int idx100 = myIdxMap[1][0][0];
        DynaComplex c001 = new DynaComplex(y[2*idx001], y[2*idx001+1]);
        DynaComplex c010 = new DynaComplex(y[2*idx010], y[2*idx010+1]);
        DynaComplex t1 = new DynaComplex(0,0);
        DynaComplex ii = new DynaComplex(0,1);
        es[0] = (0.5*n)*Math.pow(2.0, n-1)*(t1.set(c001).add(c010).getReal());
        es[1] = (0.5*n)*Math.pow(2.0, n-1)*(t1.set(c010).subtract(c001).multiply(ii).getReal());
        es[2] = (0.5*n)*Math.pow(2.0, n)*y[2*idx100];
    }
    
    @Override
    public void getSecondOrderCollectiveEvs(double[] y, DynaComplex[][] es) {
        double twonm1 = Math.pow(2.0, n-1);
        double nc2 = 0.5*n*(n-1);
        
        int idx001 = myIdxMap[0][0][1];
        int idx010 = myIdxMap[0][1][0];
        int idx100 = myIdxMap[1][0][0];
        int idx002 = myIdxMap[0][0][2];
        int idx020 = myIdxMap[0][2][0];
        int idx200 = myIdxMap[2][0][0];
        int idx110 = myIdxMap[1][1][0];
        int idx011 = myIdxMap[0][1][1];
        int idx101 = myIdxMap[1][0][1];
        
        DynaComplex c001 = new DynaComplex(y[2*idx001], y[2*idx001+1]);
        DynaComplex c010 = new DynaComplex(y[2*idx010], y[2*idx010+1]);
        DynaComplex c100 = new DynaComplex(y[2*idx100], y[2*idx100+1]);
        DynaComplex c002 = new DynaComplex(y[2*idx002], y[2*idx002+1]);
        DynaComplex c020 = new DynaComplex(y[2*idx020], y[2*idx020+1]);
        DynaComplex c200 = new DynaComplex(y[2*idx200], y[2*idx200+1]);
        DynaComplex c110 = new DynaComplex(y[2*idx110], y[2*idx110+1]);
        DynaComplex c011 = new DynaComplex(y[2*idx011], y[2*idx011+1]);
        DynaComplex c101 = new DynaComplex(y[2*idx101], y[2*idx101+1]);
        
        DynaComplex ii = new DynaComplex(0,1);
        DynaComplex t1 = new DynaComplex();
        
        t1.set(c020).add(c002).multiply(twonm1*nc2);
        es[0][0].set(c011).multiply(twonm1*n*(n-1)).add(n).add(t1).multiply(0.25);
        //es(1, 1, i) = 0.25*(twonm1*nc2*(t1.set(c020).add(c002)) ...
        //         + n*(1 + (n-1)*twonm1*c(idx011)));
        
        es[1][1].set(c011).multiply(-twonm1*n*(n-1)).add(-n).add(t1).multiply(-0.25);
        //ess(2, 2, i) = -0.25*(2^(n-1)*nc2*(c(idx020) + c(idx002)) ...
        //        - n*(1 + (n-1)*2^(n-1)*c(idx011)));
            
        
        t1.set(c020).subtract(c002).multiply(ii).multiply(nc2*0.25*twonm1);
        es[0][1].set(c100).multiply(ii).multiply(n*0.5*twonm1).add(t1);
        //ess(1, 2, i) = 1i*2^(n-3)*nc2*(c(idx020) - c(idx002)) + 1i*n*2^(n-2)*c(idx100);
   
        es[1][0].set(c100).multiply(ii).multiply(-n*0.5*twonm1).add(t1);
        //ess(2, 1, i) = 1i*2^(n-3)*nc2*(c(idx020) - c(idx002)) - 1i*n*2^(n-2)*c(idx100);
        
        t1.set(c110).add(c101).multiply(n-1);
        es[0][2].set(c010).subtract(c001).add(t1).multiply(n*0.25*twonm1);
        //ess(1, 3, i) = 2^(n-3)*n*((n-1)*(c(idx110) + c(idx101)) - c(idx001) + c(idx010));
   
        es[2][0].set(c001).subtract(c010).add(t1).multiply(n*0.25*twonm1);
        //ess(3, 1, i) = 2^(n-3)*n*((n-1)*(c(idx110) + c(idx101)) + c(idx001) - c(idx010));
        
        t1.set(c101).subtract(c110).multiply(n-1);
        es[1][2].set(c001).add(c010).multiply(-1.0).add(t1).multiply(ii).multiply(-n*0.25*twonm1);
        //ess(2, 3, i) = -1i*2^(n-3)*n*((n-1)*(c(idx101) - c(idx110)) - c(idx001) - c(idx010));
   
        es[2][1].set(c001).add(c010).add(t1).multiply(ii).multiply(-n*0.25*twonm1);
        //ess(3, 2, i) = -1i*2^(n-3)*n*((n-1)*(c(idx101) - c(idx110)) + c(idx001) + c(idx010));
        
        es[2][2].set(c200).multiply(nc2*twonm1).add(0.25*n);
        //ess(3, 3, i) = n/4 + 2^(n-1)*nc2*c(idx200);
    }

    @Override
    public void initSpinUpX(double[] y0) {
        initialize(y0, 0.5*Math.PI, 0.0, InitAngleType.CONST, InitAngleType.CONST);
    }

    @Override
    public void initialize(double[] y0, double zenith, double phase, InitAngleType zenithType,
            InitAngleType phaseType) {
        if(zenithType != InitAngleType.CONST || phaseType != InitAngleType.CONST) {
            throw new UnsupportedOperationException("SymmEval does not support the given IC");
        }
        
        int d = getDimension();
        int[] colIdx = {0,0,0};
        DynaComplex t1 = new DynaComplex(0);
        DynaComplex imag = new DynaComplex(0,1);
        for(int i = 0; i < d; ++i) {
            if(zenith == Math.PI) {
                // Spin down in z-direction
                if(colIdx[1] == 0 && colIdx[2] == 0) {
                    y0[2*i] = (1/Math.pow(2.0, n))*Math.pow(-1, colIdx[0]);
                    y0[2*i+1] = 0.0;
                } else {
                    y0[2*i] = 0.0;
                    y0[2*i+1] = 0.0;
                }
            } else if(zenith == 0) {
                // Spin up in z-direction
                if(colIdx[1] == 0 && colIdx[2] == 0) {
                    y0[2*i] = (1/Math.pow(2.0, n));
                    y0[2*i+1] = 0.0;
                } else {
                    y0[2*i] = 0.0;
                    y0[2*i+1] = 0.0;
                }
            } else if(zenith == Math.PI/2.0 && phase == 0.0) {
                // Spin up in the x-direction
                if(colIdx[0] == 0) {
                    y0[2*i] = (1/Math.pow(2.0, n));
                    y0[2*i+1] = 0.0;
                } else {
                    y0[2*i] = 0.0;
                    y0[2*i+1] = 0.0;
                }
            } else if(zenith == Math.PI/2.0 && phase == Math.PI) {
                // Spin down in the x-direction
                if(colIdx[0] == 0) {
                    y0[2*i] = (1/Math.pow(2.0, n))*Math.pow(-1, colIdx[1] + colIdx[2]);
                    y0[2*i+1] = 0.0;
                } else {
                    y0[2*i] = 0.0;
                    y0[2*i+1] = 0.0;
                }
            } else if(zenith == Math.PI/2.0 && phase == Math.PI/2) {
                // Spin up in the y-direction
                if(colIdx[0] == 0) {
                    t1.set(1,0);
                    for(int j = 0; j < colIdx[1]; ++j) {
                        t1.multiply(imag).multiply(-1);
                    }

                    for(int j = 0; j < colIdx[2]; ++j) {
                        t1.multiply(imag);
                    }
                    t1.multiply(1.0/Math.pow(2.0, n));
                    
                    y0[2*i] = t1.getReal();
                    y0[2*i+1] = t1.getImaginary();
                } else {
                    y0[2*i] = 0.0;
                    y0[2*i+1] = 0.0;
                }
            } else if(zenith == Math.PI/2.0 && phase == 3.0*Math.PI/2) {
                // Spin up in the y-direction
                if(colIdx[0] == 0) {
                    t1.set(1,0);
                    for(int j = 0; j < colIdx[1]; ++j) {
                        t1.multiply(imag);
                    }

                    for(int j = 0; j < colIdx[2]; ++j) {
                        t1.multiply(imag).multiply(-1);
                    }
                    t1.multiply(1.0/Math.pow(2.0, n));
                    
                    y0[2*i] = t1.getReal();
                    y0[2*i+1] = t1.getImaginary();
                } else {
                    y0[2*i] = 0.0;
                    y0[2*i+1] = 0.0;
                }
            } else {
                throw new UnsupportedOperationException("SymmEval does not support the given IC");
            }
            incIdx(colIdx);
        }

    }
    
    public int[][][] getIdxMap() {
        int[][][] map = new int[n+1][n+1][n+1];
        
        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < n; ++j) {
                for(int k = 0; k < n; ++k) {
                    if(i + j + k >n) {
                        map[i][j][k] = -1;
                    }
                }
            }
        }
        
        int d = getDimension();
        int[] idx = {0,0,0};
        for(int i = 0; i < d; ++i) {
            map[idx[0]][idx[1]][idx[2]] = i;
            incIdx(idx);
        }
        
        return map;
    }
    
    public void incIdx(int[] idx) {
        if(idx[0] == -1) {
            return;
        }

        ++idx[0];
        if(idx[0] + idx[1] + idx[2] > n) {
            if(idx[0] != 1) {
                idx[0] = 0;
                ++idx[1];
            } else {
                idx[0] = 0;
                idx[1] = 0;
                ++idx[2];
                if(idx[2] > n) {
                    idx[0] = -1;
                    idx[1] = -1;
                    idx[2] = -1;
                }
            }
        }
    }
}
