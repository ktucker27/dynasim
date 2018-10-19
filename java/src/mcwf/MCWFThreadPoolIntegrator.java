package mcwf;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import ode.SystemParams;
import utils.DynaComplex;

/**
 * A class for simulating a system of spin-1/2 particles using the
 * Monte Carlo wave function method. Provides a simple interface for an
 * integrator that simulates quantum trajectories in parallel and
 * organizes output
 */
public class MCWFThreadPoolIntegrator {

    private int myNumTrajectories;
    //private int myNumFinished;
    private ExecutorService myThreadPool;
    private QuantumTrajectory[] myTrajectories;
    private Runnable[] myIntegrators;
    private MCWFAggregator myAgg;
    
    public MCWFThreadPoolIntegrator(int numTrajectories, int numTimes, double timeDelta, double evDelta, SystemParams params, DynaComplex[] initialState, int numThreads, boolean debug) {
        // Parameter validation
        if(params.getN() % 2 != 0) {
            throw new UnsupportedOperationException("MCWFIntegrator requires an even number of particles");
        }
        
        if(params.getGab() != params.getGaa()) {
            throw new UnsupportedOperationException("MCWFIntegrator requires that on-site and off-site g values be equal");
        }
        
        if(initialState.length != params.getN() + 1) {
            throw new UnsupportedOperationException("Initial state has length " + initialState.length + ". Size " + (params.getN() + 1) + " expected");
        }
        
        int numEvSteps = (int)(evDelta/timeDelta);
        if(evDelta < timeDelta || Math.abs((evDelta/timeDelta) - numEvSteps) > 1.0e-10) {
            throw new UnsupportedOperationException("Expected value time delta must be an integer multiple of the time step");
        }

        myNumTrajectories = numTrajectories;
        //myNumFinished = 0;
        
        myThreadPool = Executors.newFixedThreadPool(numThreads);
        
        myTrajectories = new QuantumTrajectory[numTrajectories];
        myIntegrators = new Runnable[numTrajectories];
        
        for(int i = 0; i < numTrajectories; ++i) {
            myTrajectories[i] = new QuantumTrajectory(numTimes, numEvSteps, timeDelta, params, initialState, debug);
            //myIntegrators[i] = new MCWFFirstOrderIntegrator(myTrajectories[i]);
            myIntegrators[i] = new MCWFDelayTimeIntegrator(myTrajectories[i]);
        }
        
        myAgg = new MCWFAggregator((numTimes-1)/numEvSteps + 1, evDelta);
    }
    
    public MCWFAggregator getAggregator() {
        return myAgg;
    }
    
    /**
     * Creates trajectories and executes them by passing them to the thread pool.
     * Trajectories are stored so their results can be aggregated and output later
     */
    public void start() {
        for(int trajIdx = 0; trajIdx < myNumTrajectories; ++trajIdx) {
            myThreadPool.execute(myIntegrators[trajIdx]);
        }
    }
    
    /**
     * Blocks until all IVPs have been solved
     * 
     * @param timeout the time to wait in seconds
     * @return
     */
    public boolean waitForFinished(long timeout) {
        myThreadPool.shutdown();
        boolean success = false;
        try {
            if(myThreadPool.awaitTermination(timeout, TimeUnit.SECONDS)) {
                success = true;
                System.out.println("Processing complete");
            } else {
                myThreadPool.shutdownNow();
                if(!myThreadPool.awaitTermination(60, TimeUnit.SECONDS)) {
                    System.out.println("Thread pool failed to shut down");
                }
            }
        } catch(InterruptedException ex) {
            System.out.println("Main thread interrupted, cancelling jobs...");
            // (Re-)Cancel if current thread also interrupted
            myThreadPool.shutdownNow();
            // Preserve interrupt status
            Thread.currentThread().interrupt();
        }
        
        myAgg.aggregate(myTrajectories, 0, myNumTrajectories-1);
        
        return success;
    }
    
    public int getNumTrajectories() {
        return myNumTrajectories;
    }
    
    public QuantumTrajectory getTrajectory(int idx) {
        return myTrajectories[idx];
    }
}
