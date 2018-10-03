package mcwf;

import java.util.Random;

import ode.SystemParams;
import utils.DynaComplex;
import utils.ExpectedSpinValues;

public class QuantumTrajectory implements Runnable {
    // Simulation parameters
    private int myNumTimes;
    private int myNumEvSteps; // Number of time steps between expected value calculations
    private double myTimeDelta;
    private SystemParams myParams;
    private int myN;
    
    private double gc, gl, kl, dl, oms;
    //private double myDnj;
    
    // State indexed in [0, N+1]
    private DynaComplex[] myState, myNewState;
    private int myJ;
    
    // Output data indexed by time
    private ExpectedSpinValues[] myEvs;
    
    // Random number generator
    private Random myRng;
    
    // Temporary complex numbers for calculations
    DynaComplex t1, t2;
    
    // Constant parameters
    private final int NUM_OUTCOMES = 5;
    
    public QuantumTrajectory(int numTimes, int numEvSteps, double timeDelta, SystemParams params, DynaComplex[] initialState) {
        myNumTimes = numTimes;
        myNumEvSteps = numEvSteps;
        myTimeDelta = timeDelta;
        myParams = new SystemParams(params);
        myN = params.getN();
        
        gc = params.getGamma()*params.getFab();
        gl = params.getGamma()*(params.getFaa() - params.getFab());
        kl = 0.0; // TODO - Incoherent pumping
        dl = 0.0; // TODO - Dephasing
        oms = params.getGab()/2.0;
        //myDnj = 1.0;

        myState = new DynaComplex[params.getN() + 1];
        myNewState = new DynaComplex[params.getN() + 1];
        for(int i = 0; i < params.getN() + 1; ++i) {
            myState[i] = new DynaComplex(initialState[i]);
            myNewState[i] = new DynaComplex(0);
        }
        
        // NOTE - N is assumed to be even, and it is assumed we are starting in the highest
        //        Dicke manifold
        setJ(params.getN()/2);
        
        myEvs = new ExpectedSpinValues[(numTimes-1)/numEvSteps + 1];
        for(int i = 0; i < (numTimes-1)/numEvSteps + 1; ++i) {
            myEvs[i] = new ExpectedSpinValues();
        }
        
        myRng = new Random();
        
        t1 = new DynaComplex(0);
        t2 = new DynaComplex(0);
    }
    
    public int getNumEvTimes() {
        return myEvs.length;
    }
    
    public double getEvTimeDelta() {
        return myTimeDelta*myNumEvSteps;
    }
    
    public ExpectedSpinValues getEvs(int idx) {
        return myEvs[idx];
    }
    
    @Override
    public void run() {
        // Populate expected values at the start time
        calcEvs(myEvs[0]);
        
        double[] cdf = new double [NUM_OUTCOMES-1];
        
        for(int timeIdx = 1; timeIdx < myNumTimes; ++timeIdx) {
            // Determine the jump
            getCdf(cdf);
            int outcome = rollDice(cdf);
            
            // Update the state coefficients
            for(int i = 0; i < myNewState.length; ++i) {
                myNewState[i].set(0, 0);
            }
            
            int midx = 0;
            double ajmm = 0.0;
            double ajmp = 0.0;
            double pval = 0.0;
            switch(outcome) {
            case 0:
                // Collective decay
                for(int m = myJ; m >= -myJ + 1; --m) {
                    midx = myN/2 - m;
                    ajmm = getAjmm(myJ, m);
                    myNewState[midx + 1].set(myState[midx]).multiply(ajmm*Math.sqrt(myTimeDelta*gc));
                }
                
                break;
            case 1:
                // Individual decay, s = 0
            case 2:
                // Individual decay, s = -1
            case 3:
                // Individual decay, s = 1
                int s = 0;
                if(outcome == 2) {
                    s = -1;
                } else if(outcome == 3) {
                    s = 1;
                }
                
                for(int m = myJ; m >= -(myJ + s - 1); --m) {
                    midx = myN/2 - m;
                    ajmm = getAjmm(myJ, m);
                    if(outcome == 1) {
                        pval = getPjm0(myN, myJ, m, ajmm);
                    } else if(outcome == 2) {
                        pval = getPjmm(myN, myJ, m);
                    } else if(outcome == 3) {
                        pval = getPjmp(myN, myJ, m);
                    }
                    
                    myNewState[midx + 1].set(myState[midx]).multiply(pval*Math.sqrt(myTimeDelta*gl));
                }
                
                setJ(myJ + s);
                
                break;
            case 4:
                // No jump
                for(int m = myJ; m >= -myJ; --m) {
                    midx = myN/2 - m;
                    ajmm = getAjmm(myJ, m);
                    ajmp = getAjmp(myJ, m);
                    t1.set(oms*ajmm*ajmm,
                           -0.5*(gc*ajmm*ajmm + gl*(myN/2 + m) + kl*(myN/2 - m) + dl*myN));
                    t2.set(0, -1);
                    t1.multiply(t2).multiply(myTimeDelta).add(1.0);
                    myNewState[midx].add(t2.set(myState[midx]).multiply(t1));
                    //newc(midx,1) = newc(midx,1) + (1 - 1i*omjm*dt)*cs(midx,i);
                
                    // Coherent pumping
                    if(myParams.getOmega() > 0) {
                        if(ajmp > 0) {
                            t1.set(0, -0.5*myParams.getOmega()*ajmp*myTimeDelta);
                            myNewState[midx-1].add(t2.set(myState[midx]).multiply(t1));
                            //newc(midx-1,1) = newc(midx-1,1) - 1i*0.5*o*ajmp*dt*cs(midx,i);
                        }
                    
                        if(ajmm > 0) {
                            t1.set(0, -0.5*myParams.getOmega()*ajmm*myTimeDelta);
                            myNewState[midx+1].add(t2.set(myState[midx]).multiply(t1));
                            //newc(midx+1,1) = newc(midx+1,1) - 1i*0.5*o*ajmm*dt*cs(midx,i);
                        }
                    }
                }
                
                break;
            default:
                throw new UnsupportedOperationException("Received invalid outcome " + outcome + ". Expected in [0," + NUM_OUTCOMES + ")");
            }
            
            // Normalize the state
            double newStateNormSq = 0.0;
            for(int i = 0; i < myNewState.length; ++i) {
                newStateNormSq += Math.pow(myNewState[i].mod(), 2);
            }
            
            double newStateNorm = Math.sqrt(newStateNormSq);
            for(int i = 0; i < myState.length; ++i) {
                myState[i].set(myNewState[i]).multiply(1.0/newStateNorm);
            }
            
            // Calculate the expected values
            if(timeIdx % myNumEvSteps == 0) {
                calcEvs(myEvs[timeIdx/myNumEvSteps]);
            }
        }
    }

    private void getCdf(double[] cdf) {
        for(int i = 0; i < cdf.length; ++i) {
            cdf[i] = 0.0;
        }
        
        for(int m = myJ; m >= -myJ; --m) {
            int midx = myN/2 - m;
            double absCjmSq = Math.pow(myState[midx].mod(), 2);
            double ajmm = getAjmm(myJ, m);
            
            // Collective decay
            cdf[0] += Math.abs(ajmm*ajmm)*absCjmSq;
            
            // Individual decay
            double p0 = getPjm0(myN, myJ, m, ajmm);
            double pm = getPjmm(myN, myJ, m);
            double pp = getPjmp(myN, myJ, m);
            
            cdf[1] += Math.abs(p0*p0)*absCjmSq;
            cdf[2] += Math.abs(pm*pm)*absCjmSq;
            cdf[3] += Math.abs(pp*pp)*absCjmSq;
        }
        
        cdf[0] *= myTimeDelta*gc;
        cdf[1] *= myTimeDelta*gl;
        cdf[2] *= myTimeDelta*gl;
        cdf[3] *= myTimeDelta*gl;
        
        // Go from probabilities to CDF
        for(int i = 1; i < cdf.length; ++i) {
            cdf[i] += cdf[i-1];
        }
    }
    
    private int rollDice(double[] cdf) {
        double eps = myRng.nextDouble();
        
        int idx = 0;
        for(; idx < cdf.length; ++idx) {
            if(eps <= cdf[idx]) break;
        }
        
        return idx;
    }
    
    private void calcEvs(ExpectedSpinValues evs) {
        DynaComplex exp_p = new DynaComplex(0);
        DynaComplex exp_m = new DynaComplex(0);
        DynaComplex exp_pp = new DynaComplex(0);
        DynaComplex exp_mm = new DynaComplex(0);
        DynaComplex exp_pz = new DynaComplex(0);
        DynaComplex exp_mz = new DynaComplex(0);
        double exp_mp = 0.0;
        double exp_pm = 0.0;
        double exp_z = 0.0;
        double exp_zz = 0.0;
        for(int m = myJ; m >= -myJ; --m) {
            int midx = myN/2 - m;
            double absCjmSq = Math.pow(myState[midx].mod(), 2);
            double ajmm = getAjmm(myJ, m);
            double ajmp = getAjmp(myJ, m);
            
            if(m > -myJ) {
                exp_p.add(t1.set(myState[midx]).conjugate().multiply(myState[midx+1]).multiply(ajmm));
                exp_pz.add(t1.set(myState[midx]).conjugate().multiply(myState[midx+1]).multiply((m-1)*ajmm));
                exp_pm += absCjmSq*ajmm*ajmm;
            }
            
            if(m < myJ) {
                exp_m.add(t1.set(myState[midx]).conjugate().multiply(myState[midx-1]).multiply(ajmp));
                exp_mz.add(t1.set(myState[midx]).conjugate().multiply(myState[midx-1]).multiply((m+1)*ajmp));
                exp_mp += absCjmSq*ajmp*ajmp;
            }
            
            if(m > -myJ + 1) {
                double ajmmm = getAjmm(myJ, m-1);
                exp_pp.add(t1.set(myState[midx]).conjugate().multiply(myState[midx+2]).multiply(ajmm*ajmmm));
            }
            
            if(m < myJ - 1) {
                double ajmpp = getAjmp(myJ, m+1);
                exp_mm.add(t1.set(myState[midx]).conjugate().multiply(myState[midx-2]).multiply(ajmp*ajmpp));
            }
            
            exp_z += absCjmSq*m;
            exp_zz += absCjmSq*m*m;
        }
        
        // TODO - Double check!
        evs.setEs(0, t1.set(exp_p).add(exp_m).multiply(0.5).getReal());
        evs.setEs(1, t1.set(exp_m).subtract(exp_p).multiply(-0.5).getImaginary());
        evs.setEs(2, exp_z);
        
        evs.getEss(0, 0).set(exp_pp).add(exp_mm).add(exp_pm).add(exp_mp).multiply(0.25);
        evs.getEss(0, 1).set(exp_mm).subtract(exp_pp).add(exp_pm - exp_mp).multiply(t1.set(0, 0.25));
        evs.getEss(0, 2).set(exp_pz).add(exp_mz).multiply(0.5);
        
        evs.getEss(1, 0).set(evs.getEss(0, 1)).conjugate();
        evs.getEss(1, 1).set(exp_mm).add(exp_pp).add(-exp_mp - exp_pm).multiply(-0.25);
        evs.getEss(1, 2).set(exp_mz).subtract(exp_pz).multiply(t1.set(0, 0.5));
        
        evs.getEss(2, 0).set(evs.getEss(0, 2)).conjugate();
        evs.getEss(2, 1).set(evs.getEss(1, 2)).conjugate();
        evs.getEss(2, 2).set(exp_zz, 0);
    }
    
    private void setJ(int newj) {
        myJ = newj;
        //myDnj = getDnj(myN, myJ);
    }
    
//    private static double getDnj(int n, int j) {
//        return factorial(n)*(2*j+1)/(factorial(n/2 - j)*factorial(n/2 + j + 1));
//    }
    
    private double getAjmm(int j, int m) {
        return Math.sqrt((j+m)*(j-m+1));
    }
    
    private double getAjmp(int j, int m) {
        return Math.sqrt((j-m)*(j+m+1));
    }
    
    private double getPjm0(int n, int j, int m, double ajmm) {
        if(j == 0) {
            return 0.0;
        }
        
        return Math.sqrt((2+n)/(double)(4*j*(j+1)))*ajmm;
    }
    
    private double getPjmm(int n, int j, int m) {
        if(j == 0) {
            return 0.0;
        }
        
        return -Math.sqrt((n + 2*j + 2)*(j + m)*(j + m - 1)/(double)(4*j*(2*j + 1)));
    }
    
    private double getPjmp(int n, int j, int m) {
        return Math.sqrt((n - 2*j)*(j - m + 1)*(j - m + 2)/(double)(4*(j + 1)*(2*j + 1)));
    }
    
//    private static double factorial(int n) {
//        double ans = 1;
//        for(int i = n; i >= 1; --i) {
//            ans = ans*i;
//        }
//        return ans;
//    }
}
