package exe;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.Scanner;

import ode.CumulantAllToAllODEs;
import ode.DynaComplexODEAdapter;
import utils.DynaComplex;
import utils.SynchUtils;

public class CumulantJacobian {

    /**
     * @param args
     * @throws FileNotFoundException 
     * @throws UnsupportedEncodingException 
     */
    public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException {
        int n = 30;
        double h = 0.001;
        double gamma = 1.0;
        double delta = 0.0;
        double f = 1.0;
        double g = 10.0;
        double w = 16.75;
        
        // Get natural frequencies from Gaussian distribution
        double[] d = new double[n];
        SynchUtils.detuneGauss(delta, d);
        
        DynaComplex alpha = new DynaComplex(f, g);
        
        CumulantAllToAllODEs codes = new CumulantAllToAllODEs(n, gamma, w, 0.0, alpha, d);
        DynaComplexODEAdapter odes = new DynaComplexODEAdapter(codes);
        int dim = codes.getDimension();
        
        DynaComplex[] z0 = new DynaComplex[dim];
//        SynchUtils.initialize(z0, Math.PI/2.0, n);

        // Get initial conditions from a file
        Scanner inputStream = new Scanner(new File("/Users/kristophertucker/output/grid/glow/upper/N30/D0p0/g0p0/final_w16p75.txt"));
        inputStream.useDelimiter("\n");
        int idx = 0;
        while(inputStream.hasNext()) {
            String[] line = inputStream.next().split(",");
            z0[idx] = new DynaComplex(Double.parseDouble(line[0]), Double.parseDouble(line[1]));
            ++idx;
        }
        
        double[] y0 = new double[2*dim];
        DynaComplexODEAdapter.toReal(z0, y0);
        
        double[] yin = new double[2*dim];
        double[] yout1 = new double[2*dim];
        double[] yout2 = new double[2*dim];
        
        double[][] a = new double[2*dim][2*dim];
        
        long startTime = System.nanoTime();

        for(int i = 0; i < 2*dim; ++i) {
            yin[i] = y0[i] - h;
            odes.computeDerivatives(0.0, yin, yout1);
            
            yin[i] = y0[i] + h;
            odes.computeDerivatives(0.0, yin, yout2);

            for(int j = 0; j < 2*dim; ++j) {
                a[i][j] = (yout2[j] - yout1[j])/(2*h);
            }
            
            yin[i] = y0[i];
        }
        
        long endTime = System.nanoTime();
        
        PrintWriter writer = new PrintWriter("/Users/kristophertucker/output/temp/jacobian.txt", "UTF-8");
        for(int i = 0; i < 2*dim; ++i) {
            for(int j = 0; j < 2*dim; ++j) {
                writer.write(a[i][j] + " ");
            }
            writer.write("\n");
        }
        writer.close();

        System.out.println("Run time: " + (endTime - startTime)/1.0e9 + " seconds");
    }

}
