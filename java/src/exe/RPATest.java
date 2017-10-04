package exe;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;

import ode.CumulantAllToAllODEs;
import ode.CumulantParams;
import ode.RPAAllToAllODEs;
import utils.DynaComplex;
import utils.SynchUtils;

public class RPATest {

    /**
     * @param args
     * @throws FileNotFoundException 
     */
    public static void main(String[] args) throws FileNotFoundException {
        int n = 3;
        double gamma = 1.0;
        double delta = 5.0;
        double f = 1;
        double g = 5.0;
        
        double w = SynchUtils.getWOpt(n);
        w = 2.0;
        System.out.println("w = " + w);
        
        DynaComplex alpha = new DynaComplex(f, g);
        
        double[] d = new double[n];
        SynchUtils.detuneGauss(delta, d);
        
        CumulantParams params = new CumulantParams(n, gamma, w, delta, alpha, d);
        
        RPAAllToAllODEs odes = new RPAAllToAllODEs(params);
        CumulantAllToAllODEs codes = new CumulantAllToAllODEs(params);
        
        // Initialize everyone to spin-up along the x-direction
        int dim = odes.getDimension();
        double[] y0 = new double[dim];
        double[] yDot = new double[dim];
        double[] cyDot = new double[dim];
        for(int i = 0; i < dim; ++i) {
            y0[i] = 0.0;
        }
        
        int idx = 0;
        for(int i = 0; i < n; ++i) {
            y0[odes.getStartIdx(2) + i] = 1.0;
            for(int j = i + 1; j < n; ++j) {
                y0[odes.getStartIdx(5) + idx] = 1.0;
                ++idx;
            }
        }
                
        DynaComplex[] z0 = new DynaComplex[codes.getDimension()];
        SynchUtils.initialize(z0, Math.PI/2.0, n);
        
        // Get initial conditions from a file
        if(n > 0) {
            Scanner inputStream = new Scanner(new File("/Users/tuckerkj/output/temp/final_answer_N3_D5p0_g5p0_w3p45.txt"));
            inputStream.useDelimiter("\n");
            idx = 0;
            while(inputStream.hasNext()) {
                String[] line = inputStream.next().split(",");
                z0[idx] = new DynaComplex(Double.parseDouble(line[0]), Double.parseDouble(line[1]));
                if(z0[idx].mod() < 1.0e-10) {
                    z0[idx].set(0,0);
                }
                ++idx;
            }
            inputStream.close();
        }
        
        translate(z0, y0, codes, odes, n);
        
        for(int i = 0; i < y0.length; ++i) {
            System.out.println(y0[i]);
        }
        
        DynaComplex[] zDot = new DynaComplex[codes.getDimension()];
        for(int i = 0; i < codes.getDimension(); ++i) {
            zDot[i] = new DynaComplex(0,0);
        }
        
        odes.computeDerivatives(0, y0, yDot);
        codes.computeDerivatives(0, z0, zDot);

        translate(zDot, cyDot, codes, odes, n);
        
        for(int i = 0; i < yDot.length; ++i) {
            System.out.println(i + ": " + yDot[i] + ", " + cyDot[i] + ", " + Double.toString(yDot[i] - cyDot[i]));
        }
    }
    
    public static void translate(DynaComplex[] z, double[] y, CumulantAllToAllODEs codes, RPAAllToAllODEs odes, int n) {
        DynaComplex t1 = new DynaComplex();
        int triIdx = 0;
        int recIdx = 0;
        for(int a = 0; a < n; ++a) {
            // s^x
            y[a] = -1.0*z[codes.getStartIdx(2)/2 + a].getReal(); 
            
            // s^y
            y[odes.getStartIdx(1) + a] = 2*z[a].getImaginary();
            
            // s^z
            y[odes.getStartIdx(2) + a] = 2*z[a].getReal();
            
            for(int b = 0; b < n; ++b) {
                if(a == b) continue;
                
                if(a < b) {
                    // s^x s^x
                    y[odes.getStartIdx(3) + triIdx] = z[codes.getStartIdx(4)/2 + triIdx].getReal();
                    
                    // s^y s^y
                    y[odes.getStartIdx(4) + triIdx] = -2*t1.set(z[codes.getStartIdx(5)/2 + triIdx]).subtract(z[codes.getStartIdx(3)/2 + triIdx]).getReal();
                    
                    // s^z s^z
                    y[odes.getStartIdx(5) + triIdx] = 2*t1.set(z[codes.getStartIdx(5)/2 + triIdx]).add(z[codes.getStartIdx(3)/2 + triIdx]).getReal();

                    // s^y s^z
                    y[odes.getStartIdx(8) + recIdx] = 2*t1.set(z[codes.getStartIdx(5)/2 + triIdx]).add(z[codes.getStartIdx(3)/2 + triIdx]).getImaginary();
                    y[odes.getStartIdx(8) + b*(n-1) + a] = 2*t1.set(z[codes.getStartIdx(3)/2 + triIdx]).conjugate().add(z[codes.getStartIdx(5)/2 + triIdx]).getImaginary();

                    ++triIdx;
                }
                
                // s^x s^y
                y[odes.getStartIdx(6) + recIdx] = -2*z[codes.getStartIdx(1)/2 + recIdx].getImaginary();
                
                // s^x s^z
                y[odes.getStartIdx(7) + recIdx] = -2*z[codes.getStartIdx(1)/2 + recIdx].getReal();

                ++recIdx;
            }
        }
    }

}
