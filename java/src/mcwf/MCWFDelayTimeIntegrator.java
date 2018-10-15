package mcwf;

import org.apache.commons.math3.ode.AbstractIntegrator;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.nonstiff.ClassicalRungeKuttaIntegrator;

import ode.DynaComplexODEAdapter;

public class MCWFDelayTimeIntegrator implements Runnable {
    private QuantumTrajectory myTrajectory;
    private AbstractIntegrator myIntegrator;
    
    private final double TIME_TOL = 1.0e-10;
    
    public MCWFDelayTimeIntegrator(QuantumTrajectory trajectory) {
        myTrajectory = trajectory;
        //myIntegrator = new AdamsMoultonIntegrator(2, trajectory.getTimeDelta()*1.0e-6, trajectory.getTimeDelta(), 1.0e-3, 1.0e-5);
        myIntegrator = new ClassicalRungeKuttaIntegrator(trajectory.getTimeDelta());
        myIntegrator.addEventHandler(trajectory.getNoJumpEventHandler(), myTrajectory.getEvTimeDelta()/2, 1.0e-12, 1000);
        myIntegrator.addEventHandler(trajectory.getNoJumpStopHandler(), myTrajectory.getTimeDelta()/2, 1.0e-12, 1000);
    }
    
    @Override
    public void run() {
        int n = myTrajectory.getParams().getN();
        
        double time = 0.0;
        double finalTime = myTrajectory.getTimeDelta()*(myTrajectory.getNumTimes() - 1);
        
        double[] y0 = new double[2*(n+1)];
        double[] y = new double[2*(n+1)];
        
        double[] cdf = new double[myTrajectory.getNumOutcomes() - 1];
        
        FirstOrderDifferentialEquations odes = new DynaComplexODEAdapter(myTrajectory.getNoJumpODEs());
        
        while(time < finalTime) {
            DynaComplexODEAdapter.toReal(myTrajectory.getState(), y0);
            myIntegrator.integrate(odes, time, y0, finalTime, y);
            time = myTrajectory.getTime();
            
            // If we're not done, then it's time for a jump
            if(Math.abs(time - finalTime) <= myTrajectory.getEvTimeDelta() + TIME_TOL) {
                // We've reached the end of the integration, so calculate the final time
                myTrajectory.normalize();
                if(time < finalTime - TIME_TOL) {
                    myTrajectory.calcEvs(finalTime);
                }
                time = finalTime;
            } else if(time < finalTime) {
                myTrajectory.getCdf(cdf);
                
                // Normalize the CDF to exclude the no jump case
                double normFactor = cdf[cdf.length-1];
                for(int i = 0; i < cdf.length; ++i) {
                    cdf[i] = cdf[i]/normFactor;
                }
                
                int outcome = myTrajectory.rollDice(cdf);
                //System.out.println(time + "," + outcome);
                
                // We should never see a no jump outcome
                if(outcome == myTrajectory.getNumOutcomes() - 1) {
                    System.out.println("Time: " + time + " " + finalTime);
                    System.out.println("CDF:");
                    for(int i = 0; i < cdf.length; ++i) {
                        System.out.println(cdf[i]);
                    }
                    throw new UnsupportedOperationException("Got a jump outcome in delay time integrator");
                }
                
                myTrajectory.performJump(outcome);
            }
        }
    }
}
