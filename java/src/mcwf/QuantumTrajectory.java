package mcwf;

import java.util.Random;

import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.ode.events.EventHandler;
import org.apache.commons.math3.ode.sampling.StepHandler;
import org.apache.commons.math3.ode.sampling.StepInterpolator;

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
    
    // State indexed in [0, N+1)
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
    private final int NUM_OUTCOMES = 8;
    
    // Tolerance used to double check times
    private final double TIME_TOL = 1.0e-10;
    
    // Object to track detailed trajectory info for debug purposes
    TrajectoryStats myStats;
    
    private class NoJumpStepHandler implements StepHandler {
        private double myFinalTime;
        
        @Override
        public void init(double t0, double[] y, double tf) {
            myFinalTime = tf;
        }

        @Override
        public void handleStep(StepInterpolator interpolator, boolean isLast) throws MaxCountExceededException {
            // On the last step, set the state
            if(isLast) {
                double[] y = interpolator.getInterpolatedState();
                double t = interpolator.getInterpolatedTime();
                
                // Set and normalize the state
                DynaComplexODEAdapter.toComplex(y, myState);

                normalize();

                myTime = t;
                
                // If we're at the final time and still awaiting an expected value
                // calculation, do it now
                if(Math.abs(t - myFinalTime) < TIME_TOL && myEvIdx == myEvs.length - 1) {
                    calcEvs(t);
                }
                
                // Record the stop time
                if(myStats != null) {
                    myStats.recordStopTime(t);
                }
            }
        }

    }
    
    private class NoJumpEventHandler implements EventHandler {
        private double myTimeDelta;
        
        public NoJumpEventHandler(double timeDelta) {
            myTimeDelta = timeDelta;
        }
        
        @Override
        public void init(double t0, double[] y0, double t) {
            // Check to see if we need to compute expected values at the start
            if(Math.abs(g(t0, null)) < TIME_TOL && 
               (myEvIdx == 0 || (myEvIdx > 0 && Math.abs(t0 - myEvs[myEvIdx-1].getTime()) > TIME_TOL))) {
                DynaComplexODEAdapter.toComplex(y0, myState);
                myTime = t0;
                calcEvs(t0);
            }
        }
        
        @Override
        public EventHandler.Action eventOccurred(double t, double[] y, boolean increasing) {
            // If we are getting a redundant event at an already processed time point, continue
            if(myEvIdx > 0 && Math.abs(t - myEvs[myEvIdx-1].getTime()) < TIME_TOL) {
                return EventHandler.Action.CONTINUE;
            }
            
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
    
    public QuantumTrajectory(int numTimes, int numEvSteps, double timeDelta, SystemParams params, DynaComplex[] initialState, boolean collectStats) {
        myNumTimes = numTimes;
        myNumEvSteps = numEvSteps;
        myEvIdx = 0;
        myTimeDelta = timeDelta;
        myParams = new SystemParams(params);
        myN = params.getN();
        
        gc = params.getGamma()*params.getFab();
        gl = params.getGamma()*(params.getFaa() - params.getFab());
        kl = 0.0; // TODO - Incoherent pumping
        dl = params.getGel();
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
        
        myStats = null;
        if(collectStats) {
            myStats = new TrajectoryStats(myEvs.length);
        }
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
    
    public int getJ() {
        return myJ;
    }
    
    public DynaComplex[] getState() {
        return myState;
    }
    
    public void setState(DynaComplex[] state, int jval) {
        if(state.length != myState.length) {
            throw new UnsupportedOperationException("State provided to setState of incorrect dimension " + state.length + ". Expected: " + myState.length);
        }
        
        for(int i = 0; i < state.length; ++i) {
            myState[i].set(state[i]);
        }
        setJ(jval);
    }
    
    public void normalize() {
        double newStateNormSq = 0.0;
        for(int i = 0; i < myState.length; ++i) {
            newStateNormSq += Math.pow(myState[i].mod(), 2);
        }
        
        if(newStateNormSq == 0 || Double.isNaN(newStateNormSq)) {
            throw new UnsupportedOperationException("Found bad new state norm in normalize: " + newStateNormSq);
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
    
    public StepHandler getNoJumpStepHandler() {
        return new NoJumpStepHandler();
    }
    
    public EventHandler getNoJumpEventHandler() {
        return new NoJumpEventHandler(getEvTimeDelta());
    }
    
    public EventHandler getNoJumpStopHandler() {
        return new NoJumpStopHandler();
    }
    
    public TrajectoryStats getStats() {
        return myStats;
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
            // Dephasing, s = 0
        case 5:
            // Dephasing, s = -1
        case 6:
            // Dephasing, s = 1
            s = 0;
            if(outcome == 5) {
                s = -1;
            } else if(outcome == 6) {
                s = 1;
            }
            
            for(int m = myJ + s; m >= -(myJ + s); --m) {
                if(Math.abs(m) > myJ) continue;
                
                midx = myN/2 - m;
                if(outcome == 4) {
                    pval = getPjmz0(myN, myJ, m);
                } else if(outcome == 5) {
                    pval = getPjmzm(myN, myJ, m);
                } else if(outcome == 6) {
                    pval = getPjmzp(myN, myJ, m);
                }
                
                myNewState[midx].set(myState[midx]).multiply(pval*Math.sqrt(myTimeDelta*dl));
            }
            
            setJ(myJ + s);
            
            break;
        case 7:
            // No jump
            for(int m = myJ; m >= -myJ; --m) {
                midx = myN/2 - m;
                ajmm = getAjmm(myJ, m);
                ajmp = getAjmp(myJ, m);
                t1.set(oms*ajmm*ajmm,
                       -0.5*(gc*ajmm*ajmm + gl*(myN/2 + m) + kl*(myN/2 - m) + dl*myN/4));
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
        
        if(newStateNormSq == 0 || Double.isNaN(newStateNormSq)) {
            throw new UnsupportedOperationException("Found bad new state norm with jump: " + newStateNormSq);
        }
        
        double newStateNorm = Math.sqrt(newStateNormSq);
        for(int i = 0; i < myState.length; ++i) {
            myState[i].set(myNewState[i]).multiply(1.0/newStateNorm);
        }
        
        // Record the jump
        if(myStats != null) {
            myStats.recordJump(myTime, outcome);
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
            
            // Dephasing
            double pz0 = getPjmz0(myN, myJ, m);
            double pzm = getPjmzm(myN, myJ, m);
            double pzp = getPjmzp(myN, myJ, m);
            
            cdf[4] += Math.abs(pz0*pz0)*absCjmSq;
            cdf[5] += Math.abs(pzm*pzm)*absCjmSq;
            cdf[6] += Math.abs(pzp*pzp)*absCjmSq;
        }
        
        cdf[0] *= myTimeDelta*gc;
        cdf[1] *= myTimeDelta*gl;
        cdf[2] *= myTimeDelta*gl;
        cdf[3] *= myTimeDelta*gl;
        cdf[4] *= myTimeDelta*dl;
        cdf[5] *= myTimeDelta*dl;
        cdf[6] *= myTimeDelta*dl;
        
        if(Double.isNaN(cdf[0])) {
            throw new UnsupportedOperationException("Found NaN in CDF");
        }
        
        // Go from probabilities to CDF
        for(int i = 1; i < cdf.length; ++i) {
            cdf[i] += cdf[i-1];
            
            if(Double.isNaN(cdf[i])) {
                throw new UnsupportedOperationException("Found NaN in CDF");
            }
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
        
        // Record state information
        if(myStats != null) {
            myStats.recordState(myEvIdx - 1, time, myState, myJ);
        }
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
        return Math.sqrt((double)(j+m)*(double)(j-m+1));
    }

    private double getAjmp(int j, int m) {
        return Math.sqrt((double)(j-m)*(double)(j+m+1));
    }

    private double getPjm0(int n, int j, int m, double ajmm) {
        if(j == 0) {
            return 0.0;
        }

        return Math.sqrt((2+n)/(4*(double)j*(double)(j+1)))*ajmm;
    }

    private double getPjmm(int n, int j, int m) {
        if(j == 0) {
            return 0.0;
        }

        return -Math.sqrt((double)(n + 2*j + 2)*(double)(j + m)*(double)(j + m - 1)/(4*(double)j*(double)(2*j + 1)));
    }

    private double getPjmp(int n, int j, int m) {
        return Math.sqrt((double)(n - 2*j)*(double)(j - m + 1)*(double)(j - m + 2)/(4*(double)(j + 1)*(double)(2*j + 1)));
    }
    
    private double getPjmz0(int n, int j, int m) {
        if(j == 0) {
            return 0.0;
        }

        return Math.sqrt((2+n)/(4*(double)j*(double)(j+1)))*m;
    }

    private double getPjmzm(int n, int j, int m) {
        if(j == 0) {
            return 0.0;
        }

        return Math.sqrt((double)(n + 2*j + 2)*(double)(j - m)*(double)(j + m)/(4*(double)j*(double)(2*j + 1)));
    }

    private double getPjmzp(int n, int j, int m) {
        return Math.sqrt((double)(n - 2*j)*(double)(j - m + 1)*(double)(j + 1 + m)/(4*(double)(j + 1)*(double)(2*j + 1)));
    }
}
