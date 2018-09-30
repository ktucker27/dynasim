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
public class MCWFIntegrator {

    private int myNumTrajectories;
    private int myNumFinished;
    private double myTimeDelta;
    private int myNumTimes;
    private SystemParams myParams;
    private ExecutorService myThreadPool;
    private QuantumTrajectory[] myTrajectories;
    private MCWFAggregator myAgg;
    
    public MCWFIntegrator(int numTrajectories, int numTimes, double timeDelta, SystemParams params, DynaComplex[] initialState, int numThreads) {
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

        myNumTrajectories = numTrajectories;
        myNumFinished = 0;
        myNumTimes = numTimes;
        myTimeDelta = timeDelta;
        myParams = new SystemParams(params);
        
        myThreadPool = Executors.newFixedThreadPool(numThreads);
        
        myTrajectories = new QuantumTrajectory[numTrajectories];
        
        for(int i = 0; i < numTrajectories; ++i) {
            myTrajectories[i] = new QuantumTrajectory(numTimes, timeDelta, params, initialState);
        }
        
        myAgg = new MCWFAggregator(numTimes, timeDelta);
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
            myThreadPool.execute(myTrajectories[trajIdx]);
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
