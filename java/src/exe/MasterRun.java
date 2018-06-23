package exe;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;

import handlers.WriteHandlerMaster;

import org.apache.commons.math3.ode.nonstiff.AdamsMoultonIntegrator;

import ode.SystemParams;
import ode.DynaComplexODEAdapter;
import ode.MasterAllToAllODEs;
import utils.DynaComplex;
import utils.SynchUtils;
import utils.TPSOperator;
import utils.TPSOperator.PauliOp;

public class MasterRun {

    /**
     * @param args
     * @throws UnsupportedEncodingException 
     * @throws FileNotFoundException 
     */
    public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException {
        int n = 6;
        double h = 0.001;
        double gamma = 1.0;
        double tmax = 5.0;
        double tmin = 5.0;
        double delta = 0.0;
        double f = 1;
        double g = 0.0;
        
        DynaComplex alpha = new DynaComplex(f, g);
        
        double[] d = new double[n];
        SynchUtils.detuneGauss(delta, d);
        
        ArrayList<SystemParams> params = new ArrayList<SystemParams>();
        
        double wmin = 2.0;
        double wmax = 20.0;
        double dw = (wmax - wmin)/20;
        for(double w = wmin; w <= wmax; w += dw) {
            SystemParams p = new SystemParams(n, gamma, w, delta, alpha, d);
            params.add(p);
        }
        
        // Initialize everyone to spin-up along the x-direction
        TPSOperator rho0 = new TPSOperator(n);
        rho0.set(1.0/Math.pow(2,n));
        
        DynaComplex[] z0 = rho0.getVals();
        double[] y0 = new double[2*z0.length];
        DynaComplexODEAdapter.toReal(z0, y0);
        
        System.out.println(params.get(0).toString());
        
        long startTime = System.nanoTime();
        
        TPSOperator rho = new TPSOperator(n);
        DynaComplex t1 = new DynaComplex(0,0);
        String dir = "/Users/kristophertucker/output/master/vw/" + params.get(0).getResultsDir().getAbsolutePath() + "/";
        File fdir = new File(dir);
        fdir.mkdirs();
        PrintWriter corrWriter = new PrintWriter(dir + "corr.txt", "UTF-8");
        for(int idx = 0; idx < params.size(); ++idx) {
            SystemParams cparams = params.get(idx);
            MasterAllToAllODEs modes = new MasterAllToAllODEs(cparams);
            DynaComplexODEAdapter odes = new DynaComplexODEAdapter(modes);
            
            System.out.println("w: " + cparams.getW());
            
            //WriteHandlerMaster writeHandler = new WriteHandlerMaster("/Users/kristophertucker/output/temp/master_out.txt", n);
            AdamsMoultonIntegrator integrator = new AdamsMoultonIntegrator(2, h*1.0e-4, h, 1.0e-3, 1.0e-2);
            //integrator.addStepHandler(writeHandler);
        
            double[] y = new double[2*z0.length];
        
            integrator.integrate(odes, 0, y0, tmax, y);
            
            DynaComplexODEAdapter.toComplex(y, rho.getVals());
            
            DynaComplex sum = new DynaComplex(0,0);
            TPSOperator tmp = new TPSOperator(n);
            for(int i = 0; i < n; ++i) {
                for(int j = 0; j < n; ++j) {
                    if(i == j) continue;
                    tmp.set(rho);
                    tmp.pauliLeft(PauliOp.MINUS, j);
                    tmp.pauliLeft(PauliOp.PLUS, i);
                    tmp.trace(t1);
                    sum.add(t1);
                }
            }
            sum.multiply(1.0/(n*(n-1)));
            
            corrWriter.print(cparams.getW() + ", " + sum.getReal() + ", " + sum.getImaginary() + "\n");
            corrWriter.flush();
        }
        corrWriter.close();
        
        long endTime = System.nanoTime();
        System.out.println("Run time: " + (endTime - startTime)/1.0e9 + " seconds");
    }

}
