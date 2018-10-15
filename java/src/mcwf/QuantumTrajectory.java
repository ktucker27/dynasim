package mcwf;

import java.util.Random;

import org.apache.commons.math3.ode.events.EventHandler;

import ode.DynaComplexODEAdapter;
import ode.DynaComplexODEs;
import ode.SystemParams;
import utils.DynaComplex;
import utils.ExpectedSpinValues;

public class QuantumTrajectory {
    private int myNumTimes;
    private int myNumEvSteps;
    private int myEvIdx;
    private double myTimeDelta;
    private SystemParams myParams;
    private int myN;
    
    private double gc, gl, kl, dl, oms;
    //private double myDnj;
    
    // State indexed in [0, N+1]
    private double myTime;
    private DynaComplex[] myState, myNewState;
    private int myJ;
    
    // Output data indexed by time
    private ExpectedSpinValues[] myEvs;
    
    // No jump ODEs
    MCWFNoJumpODEs myNoJumpODEs;
    
    // Random number generator
    private Random myRng;
    
    // Temporary complex numbers for calculations
    DynaComplex t1, t2;
    
    // Constant parameters
    private final int NUM_OUTCOMES = 5;
    
    // Tolerance used to double check times
    private final double TIME_TOL = 1.0e-10;
    
    private class NoJumpEventHandler implements EventHandler {
        private double myTimeDelta;
        private double myT0;
        
        public NoJumpEventHandler(double timeDelta) {
            myTimeDelta = timeDelta;
            myT0 = 0;
        }
        
        @Override
        public void init(double t0, double[] y0, double t) {
            myT0 = t0;
            DynaComplexODEAdapter.toComplex(y0, myState);
            myTime = t0;
            if(Math.abs(g(t0, null)) < TIME_TOL) {
                calcEvs(t0);
            }
        }
        
        @Override
        public EventHandler.Action eventOccurred(double t, double[] y, boolean increasing) {
//            if(Math.abs(t - myT0) < TIME_TOL) {
//                return EventHandler.Action.CONTINUE;
//            }
            
            // Set and normalize the state
            DynaComplexODEAdapter.toComplex(y, myState);
            double newStateNormSq = 0.0;
            for(int i = 0; i < myState.length; ++i) {
                newStateNormSq += Math.pow(myState[i].mod(), 2);
            }
            
            double newStateNorm = Math.sqrt(newStateNormSq);
            for(int i = 0; i < myState.length; ++i) {
                myState[i].multiply(1.0/newStateNorm);
            }
            
            // Update the time
            myTime = t;
            
            // Calculate the expected values
            calcEvs(t);
            
            return EventHandler.Action.CONTINUE;
        }
        
        @Override
        public double g(double t, double[] y) {
            int numTimeSteps = (int)Math.floor((t+0.5*myTimeDelta)/myTimeDelta);
            double sign = 2.0*(numTimeSteps % 2) - 1.0;
            return 2.0*sign/myTimeDelta*(t - numTimeSteps*myTimeDelta);
        }
        
        @Override
        public void resetState(double t, double[] y) {
        }
    }
    
    private class NoJumpStopHandler implements EventHandler {
        private double myEps;
        
        public NoJumpStopHandler() {
            myEps = myRng.nextDouble();
        }
        
        @Override
        public void init(double t0, double[] y0, double t) {
        }
        
        @Override
        public EventHandler.Action eventOccurred(double t, double[] y, boolean increasing) {
            // Set and normalize the state
            DynaComplexODEAdapter.toComplex(y, myState);
            
            normalize();
            
            myTime = t;
            
            myEps = myRng.nextDouble();
            
            // Stop
            return EventHandler.Action.STOP;
        }
        
        @Override
        public double g(double t, double[] y) {
            DynaComplexODEAdapter.toComplex(y, myState);
            double newStateNormSq = 0.0;
            for(int i = 0; i < myState.length; ++i) {
                newStateNormSq += Math.pow(myState[i].mod(), 2);
            }
            
            return newStateNormSq - myEps;
        }
        
        @Override
        public void resetState(double t, double[] y) {
        }
    }
    
    public QuantumTrajectory(int numTimes, int numEvSteps, double timeDelta, SystemParams params, DynaComplex[] initialState) {
        myNumTimes = numTimes;
        myNumEvSteps = numEvSteps;
        myEvIdx = 0;
        myTimeDelta = timeDelta;
        myParams = new SystemParams(params);
        myN = params.getN();
        
        gc = params.getGamma()*params.getFab();
        gl = params.getGamma()*(params.getFaa() - params.getFab());
        kl = 0.0; // TODO - Incoherent pumping
        dl = 0.0; // TODO - Dephasing
        oms = params.getGab()/2.0;
        //myDnj = 1.0;

        myTime = 0.0;
        myState = new DynaComplex[params.getN() + 1];
        myNewState = new DynaComplex[params.getN() + 1];
        for(int i = 0; i < params.getN() + 1; ++i) {
            myState[i] = new DynaComplex(initialState[i]);
            myNewState[i] = new DynaComplex(0);
        }
        
        // NOTE - N is assumed to be even, and it is assumed we are starting in the highest
        //        Dicke manifold
        myNoJumpODEs = new MCWFNoJumpODEs(params);
        setJ(params.getN()/2);
        
        myEvs = new ExpectedSpinValues[(numTimes-1)/numEvSteps + 1];
        for(int i = 0; i < (numTimes-1)/numEvSteps + 1; ++i) {
            myEvs[i] = new ExpectedSpinValues();
        }
        
        myRng = new Random();
        
        t1 = new DynaComplex(0);
        t2 = new DynaComplex(0);
    }
    
    public int getNumTimes() {
        return myNumTimes;
    }
    
    public double getTimeDelta() {
        return myTimeDelta;
    }
    
    public int getNumStepsPerEv() {
        return myNumEvSteps;
    }
    
    public int getNumEvTimes() {
        return myEvs.length;
    }
    
    public double getEvTimeDelta() {
        return myTimeDelta*myNumEvSteps;
    }
    
    public SystemParams getParams() {
        return myParams;
    }
    
    public int getNumOutcomes() {
        return NUM_OUTCOMES;
    }
    
    public double getTime() {
        return myTime;
    }
    
    public DynaComplex[] getState() {
        return myState;
    }
    
    public void normalize() {
        double newStateNormSq = 0.0;
        for(int i = 0; i < myState.length; ++i) {
            newStateNormSq += Math.pow(myState[i].mod(), 2);
        }
        
        double newStateNorm = Math.sqrt(newStateNormSq);
        for(int i = 0; i < myState.length; ++i) {
            myState[i].multiply(1.0/newStateNorm);
        }
    }
    
    public ExpectedSpinValues getEvs(int idx) {
        return myEvs[idx];
    }
    
    public DynaComplexODEs getNoJumpODEs() {
        return myNoJumpODEs;
    }
    
    public EventHandler getNoJumpEventHandler() {
        return new NoJumpEventHandler(getEvTimeDelta());
    }
    
    public EventHandler getNoJumpStopHandler() {
        return new NoJumpStopHandler();
    }
    
    public void performJump(int outcome) {
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
    }
    
    public void getCdf(double[] cdf) {
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
    
    public int rollDice(double[] cdf) {
        double eps = myRng.nextDouble();
        
        int idx = 0;
        for(; idx < cdf.length; ++idx) {
            if(eps <= cdf[idx]) break;
        }
        
        return idx;
    }
    
    /*
     * Calculates and stores expected values from the current state
     */
    public void calcEvs(double time) {
        if(Math.abs(time - myEvIdx*getEvTimeDelta()) > TIME_TOL) {
            throw new UnsupportedOperationException("Time discrepancy detected in QuantumTrajectory - time: " + time + " expected: " + (myEvIdx*getEvTimeDelta()) + " traj time: " + myTime);
        }
        
        ExpectedSpinValues evs = myEvs[myEvIdx++];
        evs.setTime(time);
        
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
        myNoJumpODEs.setJ(newj);
        //myDnj = getDnj(myN, myJ);
    }
    
    //  private static double getDnj(int n, int j) {
    //  return factorial(n)*(2*j+1)/(factorial(n/2 - j)*factorial(n/2 + j + 1));
    //}

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

    //private static double factorial(int n) {
    //  double ans = 1;
    //  for(int i = n; i >= 1; --i) {
    //      ans = ans*i;
    //  }
    //  return ans;
    //}
}
