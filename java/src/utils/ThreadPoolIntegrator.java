package utils;

import integrator.IntegratorFactory;

import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;

public class ThreadPoolIntegrator {

    int myNumThreads;
    int myNumIvps;
    int myNumFinished;
    boolean myQuietMode;
    ExecutorService myThreadPool;
    ConcurrentLinkedQueue<FirstOrderIntegrator> myIntegrators;

    private class ThreadPoolIVP implements Runnable {
        FirstOrderDifferentialEquations myOdes;
        double t0;
        double[] y0;
        double t;
        double[] y;

        /**
         * @param odes the odes to integrate
         * @param t0 initial time
         * @param y0 initial solution
         * @param t final time
         * @param y solution
         */
        public ThreadPoolIVP(FirstOrderDifferentialEquations odes,
                double t0, double[] y0,
                double t, double[] y) {
            super();
            this.myOdes = odes;
            this.t0 = t0;
            this.y0 = y0;
            this.t = t;
            this.y = y;
        }

        @Override
        public void run() {
            FirstOrderIntegrator integrator = myIntegrators.remove();
            
            integrator.integrate(myOdes, t0, y0, t, y);
            
            myIntegrators.add(integrator);
            
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
    public ThreadPoolIntegrator(int numThreads, IntegratorFactory intFactory) {
        super();
        this.myNumThreads = numThreads;
        myThreadPool = Executors.newFixedThreadPool(4);
        myNumIvps = 0;
        myNumFinished = 0;
        myQuietMode = true;
        
        myIntegrators = new ConcurrentLinkedQueue<FirstOrderIntegrator>();
        for(int i = 0; i < numThreads; ++i) {
            myIntegrators.add(intFactory.newIntegrator());
        }
    }
    
    public void addIvp(FirstOrderDifferentialEquations odes, double t0, double[] y0, double t, double[] y) {
        myThreadPool.execute(new ThreadPoolIVP(odes, t0, y0, t, y));
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
