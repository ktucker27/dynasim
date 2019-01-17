package mcwf;

import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.MaxCountExceededException;

import ode.DynaComplexODEs;
import ode.SystemParams;
import utils.DynaComplex;

public class MCWFNoJumpODEs implements DynaComplexODEs {
    
    private int myN;
    private int myJ;
    private double chi;
    private double omega;
    private double gammaC;
    private double gammaS;
    private double gammaEl;
    
    private DynaComplex t1, t2;

    public MCWFNoJumpODEs(SystemParams params) {
        super();
        
        myN = params.getN();
        myJ = myN/2;
        chi = params.getGamma()*params.getGab()/2.0;
        omega = params.getOmega();
        gammaC = params.getGamma()*params.getFab();
        gammaS = params.getGamma()*(params.getFaa() - params.getFab());
        gammaEl = params.getGel();
        
        t1 = new DynaComplex(0,0);
        t2 = new DynaComplex(0,0);
    }
    
    public void setJ(int newj) {
        myJ = newj;
    }
    
    @Override
    public void computeDerivatives(double t, DynaComplex[] z, DynaComplex[] zDot)
            throws MaxCountExceededException, DimensionMismatchException {
        for(int midx = 0; midx <= myN; ++midx) {
            zDot[midx].set(0,0);
        }
        
        for(int midx = 0; midx <= myN; ++midx) {
            if(midx >= myN/2 - myJ && midx <= myN/2 + myJ) {
                int m = myN/2 - midx;
                double ajmm = getAjmm(myJ, m);
                double ajmp = getAjmp(myJ, m);
                t1.set(chi*ajmm*ajmm,
                       -0.5*(gammaC*ajmm*ajmm + gammaS*(myN/2 + m) + gammaEl*myN/4)); // + kl*(myN/2 - m) ));
                t2.set(0, -1);
                t1.multiply(t2);
                zDot[midx].add(t2.set(z[midx]).multiply(t1));
                //newc(midx,1) = newc(midx,1) + (1 - 1i*omjm*dt)*cs(midx,i);
            
                // Coherent pumping
                if(omega > 0) {
                    if(ajmp > 0) {
                        t1.set(0, -0.5*omega*ajmp);
                        zDot[midx-1].add(t2.set(z[midx]).multiply(t1));
                        //newc(midx-1,1) = newc(midx-1,1) - 1i*0.5*o*ajmp*dt*cs(midx,i);
                    }
                
                    if(ajmm > 0) {
                        t1.set(0, -0.5*omega*ajmm);
                        zDot[midx+1].add(t2.set(z[midx]).multiply(t1));
                        //newc(midx+1,1) = newc(midx+1,1) - 1i*0.5*o*ajmm*dt*cs(midx,i);
                    }
                }
            }
        }
    }

    @Override
    public int getDimension() {
        return myN + 1;
    }
    
    private double getAjmm(int j, int m) {
        return Math.sqrt((j+m)*(j-m+1));
    }
    
    private double getAjmp(int j, int m) {
        return Math.sqrt((j-m)*(j+m+1));
    }

}
