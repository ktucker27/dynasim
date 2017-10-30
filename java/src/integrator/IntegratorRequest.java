package integrator;

import java.util.ArrayList;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.events.EventHandler;
import org.apache.commons.math3.ode.sampling.StepHandler;

import utils.ODESolution;

public class IntegratorRequest {
    private FirstOrderDifferentialEquations myOdes;
    private ODESolution mySoln;
    private double[] myY0;
    private double myT0;
    private double myTF;
    
    private ArrayList<StepHandler> myStepHandlers;
    private ArrayList<EventHandler> myEventHandlers;
    
    public IntegratorRequest(FirstOrderDifferentialEquations odes, double t0, double[] y0, double tf, ODESolution soln) {
        myOdes = odes;
        mySoln = soln;
        myY0 = y0;
        myT0 = t0;
        myTF = tf;
        
        myStepHandlers = new ArrayList<StepHandler>();
        myEventHandlers = new ArrayList<EventHandler>();
    }
    
    public FirstOrderDifferentialEquations getOdes() {
        return myOdes;
    }

    public ODESolution getSoln() {
        return mySoln;
    }

    public double[] getY0() {
        return myY0;
    }

    public double getT0() {
        return myT0;
    }

    public double getTF() {
        return myTF;
    }

    public void addStepHandler(StepHandler handler) {
        myStepHandlers.add(handler);
    }
    
    public StepHandler getStepHandler(int idx) {
        return myStepHandlers.get(idx);
    }
    
    public int numStepHandlers() {
        return myStepHandlers.size();
    }
    
    public void addEventHandler(EventHandler handler) {
        myEventHandlers.add(handler);
    }
    
    public EventHandler getEventHandler(int idx) {
        return myEventHandlers.get(idx);
    }
    
    public int numEventHandlers() {
        return myEventHandlers.size();
    }
}
