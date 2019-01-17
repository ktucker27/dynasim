package ode;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;

import eval.SymmEval;
import utils.DynaComplex;

public class SymmetricODEs implements DynaComplexODEs {

    private int n;
    private double gamma;
    private double w;
    private double omega;
    private double gel;
    private double faa;
    private double fab;
    private double gaa;
    private double gab;
    
    private SymmEval myEval;
    int[][][] myIdxMap;
    
    private DynaComplex t1;
    private DynaComplex lomega, lfaa, lfab, lgaa, lgab;
    
    public SymmetricODEs(SystemParams params) {
        super();
        
        this.n= params.getN();
        this.gamma = params.getGamma();
        this.w = params.getW();
        this.omega = params.getOmega();
        this.gel = params.getGel();
        this.faa = params.getFaa();
        this.fab = params.getFab();
        this.gaa = params.getGaa();
        this.gab = params.getGab();
        
        myEval = new SymmEval(this.n);
        myIdxMap = myEval.getIdxMap();
        
        init();
    }
    
    private void init() {
        t1 = new DynaComplex(0,0);
        
        lomega = new DynaComplex(0, -0.5*omega);
        lfaa = new DynaComplex(-gamma*(faa - fab), 0);
        lfab = new DynaComplex(-gamma*fab/2.0, 0);
        lgaa = new DynaComplex(0, -gamma*gaa/2.0);
        lgab = new DynaComplex(0, -gamma*gab/2.0);
    }
    
    @Override
    public void computeDerivatives(double t, DynaComplex[] z, DynaComplex[] zDot)
            throws MaxCountExceededException, DimensionMismatchException {
        int d = getDimension();
        
        // Loop through each column of the Liouvillian adding its contribution to zDot
        for(int i = 0; i < d; ++i) {
            zDot[i].set(0,0);
        }

        int[] colIdx = {0,0,0};
        int idx = -1;
        int n1, nz, np, nm;
        for(int i = 0; i < d; ++i) {
            nz = colIdx[0];
            np = colIdx[1];
            nm = colIdx[2];
            n1 = n - nz - np - nm;
            
            // colIdx
            zDot[i].add(t1.set(z[i]).multiply(lgaa).multiply(np - nm));
            zDot[i].add(t1.set(z[i]).multiply(lfaa).multiply(nz + 0.5*(np + nm)));
            zDot[i].add(t1.set(z[i]).multiply(lfab).multiply(0.5*(nm - np)*n1 + 0.5*nz*(nm + 1 + 3*(np + 1)) 
                    + 0.5*(n1 + 1)*(np - nm) + 0.5*(nz+1)*(np + 3*nm)));
            zDot[i].add(t1.set(z[i]).multiply(-w).multiply(nz + 0.5*(np + nm)));
            zDot[i].add(t1.set(z[i]).multiply(-gel/2).multiply(np + nm));
            
            if(n1 > 0) {
                // colIdx + [1;0;0]
                idx = myIdxMap[nz+1][np][nm];
                
                zDot[idx].add(t1.set(z[i]).multiply(lgab).multiply((nz+1)*(nm-np)));
                zDot[idx].add(t1.set(z[i]).multiply(lfaa).multiply(nz+1));
                zDot[idx].add(t1.set(z[i]).multiply(lfab).multiply(0.5*(nz+1)*(nm+1+3*(np+1)) - 0.5*(nz+1)*(np - nm)));
                zDot[idx].add(t1.set(z[i]).multiply(-w).multiply(-(nz + 1)));
            }
            
            if(nz > 0) {
                // colIdx + [-1;0;0]
                idx = myIdxMap[nz-1][np][nm];
                
                zDot[idx].add(t1.set(z[i]).multiply(lgab).multiply((n1+1)*(nm-np)));
                zDot[idx].add(t1.set(z[i]).multiply(lfab).multiply(-n1*nm - n1*np - nm - np));

                // colIdx + [-1;1;0]
                idx = myIdxMap[nz-1][np+1][nm];
                
                zDot[idx].add(t1.set(z[i]).multiply(lomega).multiply(-2*(np+1)));

                // colIdx + [-1;0;1]
                idx = myIdxMap[nz-1][np][nm+1];
                
                zDot[idx].add(t1.set(z[i]).multiply(lomega).multiply(2*(nm+1)));
            }
            
            if(nz > 0 && n1 > 0) {
                // colIdx + [-1;1;1]
                idx = myIdxMap[nz-1][np+1][nm+1];
                
                zDot[idx].add(t1.set(z[i]).multiply(lfab).multiply(-4*(np+1)*(nm+1)));
            }
            
            if(nz > 1) {
                // colIdx + [-2;1;1]
                idx = myIdxMap[nz-2][np+1][nm+1];
                
                zDot[idx].add(t1.set(z[i]).multiply(lfab).multiply(-4*(np+1)*(nm+1)));
            }
            
            if(np > 0 && nm > 0) {
                // colIdx + [1;-1;-1]
                idx = myIdxMap[nz+1][np-1][nm-1];
                
                zDot[idx].add(t1.set(z[i]).multiply(lfab).multiply((n1+1)*(nz+1)));
            
                // colIdx + [2;-1;-1]
                idx = myIdxMap[nz+2][np-1][nm-1];
                
                zDot[idx].add(t1.set(z[i]).multiply(lfab).multiply(-(nz+1)*(nz+2)));
            }
            
            if(nm > 0) {
                // colIdx + [1;0;-1]
                idx = myIdxMap[nz+1][np][nm-1];
                
                zDot[idx].add(t1.set(z[i]).multiply(lomega).multiply(nz+1));
            }
            
            if(np > 0) {
                // colIdx + [1;-1;0]
                idx = myIdxMap[nz+1][np-1][nm];
                
                zDot[idx].add(t1.set(z[i]).multiply(lomega).multiply(-(nz+1)));
            }
            
            myEval.incIdx(colIdx);
        }
    }

    @Override
    public int getDimension() {
        return myEval.getDimension();
    }

}
