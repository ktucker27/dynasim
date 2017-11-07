package integrator;


import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Semaphore;
import java.util.concurrent.TimeUnit;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;

import utils.ODESolution;

public class ThreadPoolIntegrator {

    private int myNumIvps;
    private int myNumFinished;
    private boolean myQuietMode;
    private ExecutorService myThreadPool;
    private ConcurrentLinkedQueue<FirstOrderIntegrator> myIntegrators;
    private Semaphore myMutex;
    private IntegratorReducer myReducer;

    private class ThreadPoolIVP implements Runnable {
        private IntegratorRequest myRequest;
        
        /**
         * @param odes the odes to integrate
         * @param t0 initial time
         * @param y0 initial solution
         * @param t final time
         * @param soln solution
         */
        public ThreadPoolIVP(IntegratorRequest request) {
            super();
            myRequest = request;
        }

        @Override
        public void run() {
            // Grab an integrator off the queue
            FirstOrderIntegrator integrator = myIntegrators.remove();
            
            double[] y = new double[myRequest.getOdes().getDimension()];
            try {
                // Setup step handlers
                integrator.clearStepHandlers();
                for(int i = 0; i < myRequest.numStepHandlers(); ++i) {
                    integrator.addStepHandler(myRequest.getStepHandler(i));
                }
                
                // Setup event handlers
                integrator.clearEventHandlers();
                for(int i = 0; i < myRequest.numEventHandlers(); ++i) {
                    integrator.addEventHandler(myRequest.getEventHandler(i), Double.POSITIVE_INFINITY, 1.0e-12, 100);
                }
                
                // Perform integration
                integrator.integrate(myRequest.getOdes(), 
                                     myRequest.getT0(),
                                     myRequest.getY0(),
                                     myRequest.getTF(),
                                     y);
                
                // Set solution
                myRequest.getSoln().setSolution(y);
            } catch(Exception ex) {
                System.err.println(ex.getMessage());
                ex.printStackTrace();
            } 
            
            // Perform reduction step
            if(myReducer != null) {
                try {
                    myMutex.acquire();
                    myReducer.reduce(myRequest.getSoln());
                } catch (InterruptedException ex) {
                    System.err.println(ex.getMessage());
                    ex.printStackTrace();
                } finally {
                    myMutex.release();
                }
            }
            
            // Return the integrator to the queue
            myIntegrators.add(integrator);
            
            // Post some status
            if(!myQuietMode) {
                int numFinished = finished();
                if(numFinished % 100 == 0) {
                    System.out.println("Finished " + numFinished + " out of " + getNumTodo());
                }
            }
        }
    }
        
    /**
     * @param numThreads the number of threads to use
     * @param intFactory factory to generate the integrators that will be used
     */
    public ThreadPoolIntegrator(int numThreads, IntegratorFactory intFactory, IntegratorReducer reducer) {
        super();
        myThreadPool = Executors.newFixedThreadPool(numThreads);
        myNumIvps = 0;
        myNumFinished = 0;
        myQuietMode = true;
        myReducer = reducer;
        
        myIntegrators = new ConcurrentLinkedQueue<FirstOrderIntegrator>();
        for(int i = 0; i < numThreads; ++i) {
            myIntegrators.add(intFactory.newIntegrator());
        }
        
        myMutex = new Semaphore(1);
    }
    
    public void addIvp(FirstOrderDifferentialEquations odes, double t0, double[] y0, double t, ODESolution soln) {
        IntegratorRequest request = new IntegratorRequest(odes, t0, y0, t, soln);
        myThreadPool.execute(new ThreadPoolIVP(request));
        todo();
    }
    
    public void addIvp(IntegratorRequest request) {
        myThreadPool.execute(new ThreadPoolIVP(request));
        todo();
    }
    
    public void setQuietMode(boolean quietMode) {
        myQuietMode = quietMode;
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
        
        return success;
    }
    
    private synchronized int finished() {
        return ++myNumFinished;
    }
    
    private synchronized void todo() {
        ++myNumIvps;
    }
    
    private synchronized int getNumTodo() {
        return myNumIvps;
    }
}
